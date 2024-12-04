import os
import sys

# Remove current directory from path to force using installed package
sys.path = [p for p in sys.path if not os.path.abspath(".") == p]
