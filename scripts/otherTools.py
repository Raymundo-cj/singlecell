"""
@Date: 2020-06-05 22:08:50
LastEditors: liuzj
LastEditTime: 2021-01-29 13:20:18
@Description: 无法归类的工具
@Author: liuzj
FilePath: /jpy_tools/otherTools.py
"""
import os
from functools import partial as _partial
import sh
import pandas as pd
import numpy as np
from loguru import logger
from io import StringIO
import sys
from tempfile import NamedTemporaryFile
from threading import Thread
import matplotlib.pyplot as plt
import matplotlib as mpl
import patchworklib as pw
import seaborn as sns
from matplotlib.widgets import PolygonSelector
from matplotlib.path import Path

from typing import (
    List,
    Optional,
    Union,
    Sequence,
    Literal,
    Any,
    Tuple,
    Iterator,
    Mapping,
    Callable,
    Dict,
)


class F(_partial):
    """
    Python Pipe. e.g.`range(10) | F(filter, lambda x: x % 2) | F(sum)`
    """

    def __call__(self, *args, **keywords):
        args_iter = iter(args)
        return self.func(
            *map(lambda arg: (next(args_iter) if arg == ... else arg), self.args),
            *args_iter,
            **{**self.keywords, **keywords},
        )

    def __ror__(self, other):
        return self(other)

    def __rrshift__(self, other):
        return self(other)


def setSeed(seed=0):
    import os
    import numpy as np
    import random
    import rpy2.robjects as ro
    import torch
    
    R = ro.r

    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)
    R("set.seed")(seed)
    
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.enabled = False


class Capturing(list):
    "Capture std output"

    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stdout = self._stdout


def mkdir(dirPath):
    try:
        sh.mkdir(dirPath)
    except:
        logger.warning(f"{dirPath} existed!!")
