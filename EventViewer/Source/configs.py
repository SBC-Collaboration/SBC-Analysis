

import logging
from PIL import PngImagePlugin
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

from .PopUpHandler import PopUpHandler


PngImagePlugin.MAX_TEXT_CHUNK = 2000
# verbosity = logging.DEBUG
verbosity = logging.INFO

DEFAULT_FONT=("Arial", 20)


logger = logging.getLogger('ped')
formatter = logging.Formatter(fmt='%(asctime)s [%(levelname)s] %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
console = logging.StreamHandler()
console.setFormatter(formatter)

logger.addHandler(console)
logger.setLevel(verbosity)

pop_up = PopUpHandler()
pop_up.setLevel(logging.ERROR)
logger.addHandler(pop_up)