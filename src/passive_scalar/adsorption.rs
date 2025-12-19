use super::Node;
use crate::prelude_crate::*;

pub enum AdsorptionModel {
    Langmuir {
        adsorption_rate: Float,
        desorption_rate: Float,
        q_max: Float,
    },
}

pub struct Parameters {
    adsorption_model: AdsorptionModel,
    initial_q_ads: Float,
}

impl Parameters {
    pub fn new(adsorption_model: AdsorptionModel, initial_q_ads: Float) -> Self {
        Parameters {
            adsorption_model,
            initial_q_ads,
        }
    }
}

impl Parameters {
    fn get_adsorption_model(&self) -> &AdsorptionModel {
        &self.adsorption_model
    }

    pub(super) fn get_initial_q_ads(&self) -> Float {
        self.initial_q_ads
    }
}

impl Node {
    fn get_q_ads(&self) -> Option<Float> {
        *self.q_ads.read().unwrap()
    }

    fn set_q_ads(&self, q_ads: Float) {
        let mut q_ads_guard = self.q_ads.write().unwrap();
        *q_ads_guard = Some(q_ads);
    }

    fn get_adsorption_parameters(&self) -> &Option<Parameters> {
        &self.adsorption_parameters
    }

    pub(super) fn compute_adsorption_source_value(&self) -> Option<Float> {
        if let Some(ads_params) = self.get_adsorption_parameters() {
            if self.is_bounce_back_node() {
                let scalar_value = self.get_scalar_value();
                let q_ads = self.get_q_ads().unwrap();
                let r = match ads_params.get_adsorption_model() {
                    AdsorptionModel::Langmuir {
                        adsorption_rate,
                        desorption_rate,
                        q_max,
                    } => adsorption_rate * scalar_value * (q_max - q_ads) - desorption_rate * q_ads,
                };
                let source_value = -r;
                self.update_q_ads(r);
                Some(source_value)
            } else {
                Some(0.0)
            }
        } else {
            None
        }
    }

    fn update_q_ads(&self, r: Float) {
        if self.is_bounce_back_node() {
            let q_ads_old = self.get_q_ads().unwrap_or(0.0);
            let mut q_ads_new = q_ads_old + r * DELTA_T;
            let q_max = match self
                .get_adsorption_parameters()
                .as_ref()
                .unwrap()
                .get_adsorption_model()
            {
                AdsorptionModel::Langmuir { q_max, .. } => *q_max,
            };
            if q_ads_new < 0.0 {
                q_ads_new = 0.0;
            } else if q_ads_new > q_max {
                q_ads_new = q_max;
            }
            self.set_q_ads(q_ads_new);
        }
    }
}
