import os
import sys
import time
import logging

from pylenm2.utils import constants as c


# def log_separator(logger, msg=f"\n\n\n{'-'*50}\n\n\n", *args, **kwargs):
#         # Get filehandler
#         original_fh_formatters = [(i, handler.formatter) for i, handler in enumerate(logger.handlers) if isinstance(handler, logging.FileHandler)]
             
#         # Log separators without any extra formatting other than the message
#         for i, orig_formatter in original_fh_formatters:
#              logger.handlers[i].setFormatter(logging.Formatter('%(message)s'))
#              logger.info(msg, *args, **kwargs)
#              logger.handlers[i].setFormatter(orig_formatter)



def setup_logging(
        module_name,
        level=None,
        logfile_dir=c.LOGFILE_DIR,
        get_stream_logs=True,
        main_module=False,
    ):

    # Make log directory if it doesn't exist
    os.makedirs(logfile_dir, exist_ok=True)

    # define default logging levels
    if level is None:
        default_level = logging.INFO
        filehandler_level = logging.DEBUG
        streamhandler_level = logging.ERROR if get_stream_logs else None
    else:
        default_level = level
        filehandler_level = level
        streamhandler_level = level if get_stream_logs else None

    # Create a logger object
    logger = logging.getLogger(module_name)
    logger.setLevel(default_level)  # logging.INFO by default

    if not logger.handlers:  # Avoid adding multiple handlers if already exist
        # Create file handler which logs messages in a file. One log file is used for each day.
        current_date = time.strftime("%Y-%m-%d", time.localtime())
        fh = logging.FileHandler(f"{logfile_dir}/{current_date}.log")
        fh.setLevel(filehandler_level)    # logging.DEBUG by default

        # Create console handler with a higher log level
        if get_stream_logs:
            ch = logging.StreamHandler(sys.stdout)
            ch.setLevel(streamhandler_level)    # logging.ERROR by default

        # Create formatter and add it to the handlers
        formatter = logging.Formatter(c.LOG_FORMAT)
        fh.setFormatter(formatter)
        if get_stream_logs:
            ch.setFormatter(formatter)

        # Add the handlers to the logger
        logger.addHandler(fh)
        if get_stream_logs:
            logger.addHandler(ch)

    if main_module:
        logger.info(c.LOG_SEPARATOR)

    return logger