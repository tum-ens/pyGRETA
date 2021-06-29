import logging
import sys

logger = logging.getLogger()
# logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s;%(levelno)s %(processName)s; %(filename)s: %(funcName)s - %(levelname)s: %(message)s')
#formatter = logging.Formatter('%Y-%m-%d %H:%M:%S %(message)s')

file_handler = logging.FileHandler('output.log')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)    # Logging to file

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)    # Logging to console
# logger.addHandler(logging.StreamHandler(sys.stderr))  # Logging to console as error in red
