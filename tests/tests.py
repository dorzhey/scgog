import sys
from scgog.module1 import add_one

import pytest

def test_module1():
    assert add_one(3) == 4