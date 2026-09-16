"""Enable ``python -m zodiac_dating``."""
import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
