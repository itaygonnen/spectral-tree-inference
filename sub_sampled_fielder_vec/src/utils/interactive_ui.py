"""Interactive UI components for STDR launcher.

Provides colored text, ASCII art logo, and input helpers for the interactive launcher.
"""
import os
import sys
from typing import Optional, List, Dict, Any, Sequence, Tuple


# ANSI color codes for gradient effects
class Colors:
    """ANSI escape codes for colored terminal output."""
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    MAGENTA = '\033[95m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    RESET = '\033[0m'

    # 24-bit colour where the terminal advertises it, otherwise the 6x6x6 cube of
    # the 256-colour palette — which Terminal.app and every other ANSI terminal
    # render correctly. The eight named codes above are too few for a gradient.
    TRUECOLOR = os.environ.get("COLORTERM", "").lower() in ("truecolor", "24bit")

    @staticmethod
    def rgb(r: int, g: int, b: int) -> str:
        """Foreground escape for one RGB triple, at the best depth available."""
        r, g, b = (max(0, min(255, int(v))) for v in (r, g, b))
        if Colors.TRUECOLOR:
            return f"\033[38;2;{r};{g};{b}m"
        cube = 16 + 36 * round(r / 51) + 6 * round(g / 51) + round(b / 51)
        return f"\033[38;5;{cube}m"

    @staticmethod
    def ramp(stops: Sequence[Tuple[int, int, int]], n: int) -> List[str]:
        """``n`` escapes interpolated along the piecewise-linear path through ``stops``."""
        if n <= 0:
            return []
        if n == 1 or len(stops) == 1:
            return [Colors.rgb(*stops[0])]
        out = []
        span = len(stops) - 1
        for i in range(n):
            pos = i / (n - 1) * span
            k = min(int(pos), span - 1)
            t = pos - k
            a, b = stops[k], stops[k + 1]
            out.append(Colors.rgb(*(a[j] + (b[j] - a[j]) * t for j in range(3))))
        return out

    @staticmethod
    def gradient(text: str, colors: List[str]) -> str:
        """Apply color gradient to text line-by-line."""
        lines = text.split('\n')
        if len(lines) <= 1:
            return colors[0] + text + Colors.RESET

        colored_lines = []
        for i, line in enumerate(lines):
            color_idx = int((i / (len(lines) - 1)) * (len(colors) - 1))
            colored_lines.append(colors[color_idx] + line + Colors.RESET)

        return '\n'.join(colored_lines)


# One glyph per letter, six rows each, padded to a fixed width so the letters can
# be spaced out and coloured independently.
_STDR_GLYPHS = {
    "S": ["███████╗",
          "██╔════╝",
          "███████╗",
          "╚════██║",
          "███████║",
          "╚══════╝"],
    "T": ["████████╗",
          "╚══██╔══╝",
          "   ██║   ",
          "   ██║   ",
          "   ██║   ",
          "   ╚═╝   "],
    "D": ["██████╗ ",
          "██╔══██╗",
          "██║  ██║",
          "██║  ██║",
          "██████╔╝",
          "╚═════╝ "],
    "R": ["██████╗ ",
          "██╔══██╗",
          "██████╔╝",
          "██╔══██╗",
          "██║  ██║",
          "╚═╝  ╚═╝"],
}

# Gradient stops, walked **top to bottom** across the glyph rows. Six stops
# interpolated over six rows means every row gets its own colour rather than three
# banded pairs.
_LOGO_STOPS = (
    (34, 211, 238),    # cyan
    (56, 152, 245),    # sky
    (99, 102, 241),    # indigo
    (147, 92, 246),    # violet
    (205, 74, 235),    # fuchsia
    (244, 114, 182),   # pink
)

# Gap between glyphs, wide enough that the letters read as "S T D R".
_LOGO_GAP = "    "


def _letter_spaced(word: str) -> str:
    """``'sub-sampled'`` -> ``'s u b - s a m p l e d'``."""
    return " ".join(word)


def print_logo():
    """Letter-spaced 'sub-sampled' over STDR, shaded down the same gradient.

    The block letters shade **vertically** (one colour per glyph row) and the
    headline shades horizontally along the identical ramp — on a single line
    there is no vertical axis to use, and sharing the ramp ties the two together.
    """
    pad = "    "
    rows = _STDR_GLYPHS["S"]

    head = _letter_spaced("sub-sampled")
    head_cols = Colors.ramp(_LOGO_STOPS, len(head))
    print()
    print(pad + Colors.BOLD
          + "".join(c + ch for c, ch in zip(head_cols, head)) + Colors.RESET)

    for row, col in enumerate(Colors.ramp(_LOGO_STOPS, len(rows))):
        line = _LOGO_GAP.join(_STDR_GLYPHS[ch][row] for ch in "STDR")
        print(f"{pad}{col}{line}{Colors.RESET}")

    print(f"\n{pad}{Colors.rgb(*_LOGO_STOPS[2])}Spectral Tree Recovery"
          f"{Colors.RESET}")
    print()


def print_header(text: str):
    """Print section header."""
    print(f"\n{Colors.BOLD}{Colors.YELLOW}{text}{Colors.RESET}")
    print("─" * len(text))


def print_option(key: str, description: str, highlight: bool = False):
    """Print menu option."""
    if highlight:
        print(f"  {Colors.GREEN}[{key}]{Colors.RESET} {Colors.BOLD}{description}{Colors.RESET}")
    else:
        print(f"  {Colors.CYAN}[{key}]{Colors.RESET} {description}")


def print_cached_matrix(index: int, metadata: Dict[str, Any]):
    """Print cached matrix info in compact format."""
    n = metadata.get('n_taxa', '?')
    L = metadata.get('seq_len', '?')
    mu = metadata.get('mutation_rate', '?')
    model = metadata.get('tree_model', '?')

    print(f"  {Colors.CYAN}[{index}]{Colors.RESET} n={n}, L={L}, μ={mu}, {model}")


def print_config_summary(config: Dict[str, Any], prefix: str = "  →"):
    """Print comprehensive config summary showing matrix and experiment parameters."""
    # Matrix parameters
    n = config.get('n_taxa', config.get('taxa_values', [None])[0])
    L = config.get('seq_len', config.get('sequence_length_values', [None])[0])
    mu = config.get('mutation_rate', '?')
    model = config.get('tree_model', '?')

    # Experiment parameters
    bootstrap_reps = config.get('bootstrap_reps', '?')
    sampling_method = config.get('sampling_method', '?')
    display_mode = config.get('display_mode', '?')

    # P-values info
    p_values = config.get('p_values', [])
    if p_values:
        p_min = min(p_values)
        p_max = max(p_values)
        p_count = len(p_values)
        p_info = f"{p_count} points [{p_min:.1e}-{p_max:.1e}]"
    else:
        p_info = "?"

    # Matrix line
    print(f"{prefix} {Colors.BOLD}Matrix:{Colors.RESET} n={n}, L={L}, μ={mu}, {model}")

    # Experiment line
    print(f"{prefix} {Colors.BOLD}Experiment:{Colors.RESET} {bootstrap_reps} bootstraps, p={p_info}, {display_mode} mode")

    # Sampling line
    if sampling_method == 'leveraged':
        theta = config.get('sampling_theta', '?')
        rank = config.get('sampling_target_rank', '?')
        print(f"{prefix} {Colors.BOLD}Sampling:{Colors.RESET} leveraged (θ={theta}, rank={rank})")
    else:
        print(f"{prefix} {Colors.BOLD}Sampling:{Colors.RESET} {sampling_method}")


def get_input(prompt: str, default: Optional[str] = None, color: str = Colors.CYAN) -> str:
    """Get user input with optional default value."""
    if default is not None:
        full_prompt = f"{prompt} {Colors.BOLD}[{default}]{Colors.RESET}: "
    else:
        full_prompt = f"{prompt}: "

    while True:
        try:
            user_input = input(color + full_prompt + Colors.RESET).strip()
            return user_input if user_input else default
        except (KeyboardInterrupt, EOFError):
            print(f"\n{Colors.RED}Interrupted by user{Colors.RESET}")
            sys.exit(0)
        except UnicodeDecodeError:
            # e.g. a key pressed while the keyboard is on a non-Latin layout: the byte
            # never reaches us as text. Re-ask instead of dying mid-run.
            print(f"{Colors.YELLOW}Could not read that keystroke (non-UTF-8 input) -- "
                  f"check the keyboard layout and type it again{Colors.RESET}")


def get_choice(prompt: str, valid_choices: List[str], case_sensitive: bool = False) -> str:
    """Get user choice from valid options."""
    while True:
        choice = get_input(prompt)

        if choice is None:
            continue

        if not case_sensitive:
            choice = choice.lower()
            valid_choices = [c.lower() for c in valid_choices]

        if choice in valid_choices:
            return choice

        print(f"{Colors.RED}Not an option, try again. Valid choices: {', '.join(valid_choices)}{Colors.RESET}")


def get_multi_choice(prompt: str, valid_choices: List[str], case_sensitive: bool = False) -> List[str]:
    """
    Get multiple user choices from valid options (comma-separated).

    Examples:
        >>> get_multi_choice("Choice", ["1", "2", "3", "n", "q"])
        User enters: "1,3"
        Returns: ["1", "3"]

        >>> get_multi_choice("Choice", ["1", "2", "3", "n", "q"])
        User enters: "n"
        Returns: ["n"]

    Args:
        prompt: Question to ask user
        valid_choices: List of valid options
        case_sensitive: Whether choices are case-sensitive

    Returns:
        List of selected choices (single item if no comma, multiple if comma-separated)
    """
    while True:
        choice_str = get_input(prompt)

        if choice_str is None:
            continue

        # Split by comma and strip whitespace
        choices = [c.strip() for c in choice_str.split(',')]

        if not case_sensitive:
            choices = [c.lower() for c in choices]
            valid_choices_lower = [c.lower() for c in valid_choices]
        else:
            valid_choices_lower = valid_choices

        # Validate all choices
        invalid = [c for c in choices if c not in valid_choices_lower]
        if invalid:
            print(f"{Colors.RED}Invalid choices: {', '.join(invalid)}. Valid: {', '.join(valid_choices)}{Colors.RESET}")
            continue

        return choices


def get_menu_choice(prompt: str, options: List[str], default_index: int = 0) -> str:
    """
    Present numbered menu and get user choice.

    Args:
        prompt: Question to ask user
        options: List of option strings to choose from
        default_index: Index of default option (0-based)

    Returns:
        Selected option string
    """
    # Display menu
    print(f"\n{Colors.CYAN}{prompt}{Colors.RESET}")
    for i, option in enumerate(options, 1):
        if i - 1 == default_index:
            print(f"  {Colors.GREEN}[{i}]{Colors.RESET} {Colors.BOLD}{option}{Colors.RESET} (default)")
        else:
            print(f"  {Colors.CYAN}[{i}]{Colors.RESET} {option}")

    # Get choice
    while True:
        choice = get_input("Choice", default=str(default_index + 1))

        if choice is None:
            return options[default_index]

        try:
            choice_idx = int(choice) - 1
            if 0 <= choice_idx < len(options):
                return options[choice_idx]
            else:
                print(f"{Colors.RED}Invalid choice. Enter 1-{len(options)}{Colors.RESET}")
        except ValueError:
            print(f"{Colors.RED}Invalid input. Enter a number 1-{len(options)}{Colors.RESET}")


def confirm(prompt: str = "Continue?", default: bool = True) -> bool:
    """Ask for yes/no confirmation."""
    default_str = "Y/n" if default else "y/N"

    while True:
        response = get_input(f"{prompt} [{default_str}]", default="y" if default else "n")

        if response is None:
            return default

        response_lower = response.lower()
        if response_lower in ['y', 'yes']:
            return True
        elif response_lower in ['n', 'no']:
            return False
        else:
            print(f"{Colors.RED}Not an option, try again. Enter 'y' or 'n'{Colors.RESET}")


def print_error(message: str):
    """Print error message."""
    print(f"{Colors.RED}✗ Error: {message}{Colors.RESET}")


def print_success(message: str):
    """Print success message."""
    print(f"{Colors.GREEN}✓ {message}{Colors.RESET}")


def print_warning(message: str):
    """Print warning message."""
    print(f"{Colors.YELLOW}⚠ {message}{Colors.RESET}")


def clear_screen():
    """Clear terminal screen (optional, for cleaner UX)."""
    import os
    os.system('cls' if os.name == 'nt' else 'clear')


def print_divider():
    """Print a visual divider."""
    print(f"{Colors.CYAN}{'═' * 60}{Colors.RESET}")
