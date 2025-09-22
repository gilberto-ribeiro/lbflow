from enum import Enum
from dataclasses import dataclass
from typing import List, Tuple, Dict
from PIL import Image


class Color(Enum):
    WHITE = 0
    BLACK = 1
    RED = 2
    GREEN = 3
    BLUE = 4
    YELLOW = 5
    CYAN = 6
    MAGENTA = 7


@dataclass
class InitialScalar:
    file_name: str
    header: str
    initial_values: Dict[Color, float]
    default_value: float = 0.0


class DataImage:

    def __init__(self, file_name: str):
        self.image = Image.open(file_name).convert('RGB')
        self.width, self.height = self.image.size
        self.colors: List[Tuple[int, int, Color]] = self._set_colors()

    def _set_colors(self) -> List[Tuple[int, int, Color]]:
        colors: List[Tuple[int, int, Color]] = []
        for y in range(self.height):
            for x in range(self.width):
                r, g, b = self.image.getpixel((x, y))
                color: Color = self.convert_rgb_to_color(r, g, b)
                colors.append((x, y, color))
        return colors

    def write_bounce_back_map_file(self) -> None:
        limit: int = self.width - 1
        with open("map.xyz", "w") as f:
            for x, _, color in self.colors:
                code: str
                match color:
                    case Color.BLACK:
                        code = "1"
                    case _:
                        code = "0"
                f.write(code)
                if x == limit:
                    f.write("\n")
                else:
                    f.write(" ")

    def write_initial_scalar_file(self, initial_scalar: InitialScalar) -> None:
        with open(initial_scalar.file_name, "w") as f:
            f.write(initial_scalar.header + "\n")
            for _, _, color in self.colors:
                value: float = initial_scalar.initial_values.get(
                    color, initial_scalar.default_value)
                f.write(f"{value}\n")

    @staticmethod
    def convert_rgb_to_color(r: int, g: int, b: int) -> Color:
        r_bool: bool = r > 127
        g_bool: bool = g > 127
        b_bool: bool = b > 127
        match (r_bool, g_bool, b_bool):
            case (True, True, True):
                return Color.WHITE
            case (False, False, False):
                return Color.BLACK
            case (True, False, False):
                return Color.RED
            case (False, True, False):
                return Color.GREEN
            case (False, False, True):
                return Color.BLUE
            case (True, True, False):
                return Color.YELLOW
            case (False, True, True):
                return Color.CYAN
            case (True, False, True):
                return Color.MAGENTA
            case _:
                raise ValueError("Invalid RGB values")


def main() -> None:
    pass


if __name__ == "__main__":
    main()
