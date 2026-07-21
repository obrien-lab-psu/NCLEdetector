"""
Internal logging helper for NCLEdetector.

Usage within a class ``__init__``::

    from NCLEdetector._logging import setup_logger
    self.logger = setup_logger('ClassName', outdir=self.outdir, ID=self.ID, log_level=log_level)

Users of the package can control verbosity at the top-level ``NCLEdetector`` logger::

    import logging
    logging.getLogger('NCLEdetector').setLevel(logging.WARNING)   # suppress INFO
    logging.getLogger('NCLEdetector').setLevel(logging.DEBUG)     # enable DEBUG

"""
import logging
import os


def setup_logger(
    name: str,
    outdir: str = None,
    ID: str = '',
    log_level: int = logging.INFO,
) -> logging.Logger:
    """
    Create and return a named logger for an NCLEdetector class.

    Parameters
    ----------
    name : str
        Short class name used as the logger suffix (e.g. ``'GaussianEntanglement'``).
    outdir : str, optional
        If provided, a ``FileHandler`` writing to ``<outdir>/<ID>.log`` (or
        ``<outdir>/<name>.log`` when *ID* is empty) is added alongside the
        stream handler.
    ID : str, optional
        Identifier used as the log file stem when *outdir* is given.
    log_level : int, optional
        Logging level for this logger, by default ``logging.INFO``.

    Returns
    -------
    logging.Logger
        Configured logger under the ``NCLEdetector.<name>`` hierarchy.
    """
    logger = logging.getLogger(f'NCLEdetector.{name}')
    logger.setLevel(log_level)

    # Avoid adding duplicate handlers when the class is re-instantiated
    if not logger.handlers:
        fmt = logging.Formatter(
            '%(asctime)s [%(name)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
        )

        sh = logging.StreamHandler()
        sh.setFormatter(fmt)
        logger.addHandler(sh)

        if outdir:
            os.makedirs(outdir, exist_ok=True)
            logfilename = f'{ID}.log' if ID else f'{name}.log'
            fh = logging.FileHandler(os.path.join(outdir, logfilename))
            fh.setFormatter(fmt)
            logger.addHandler(fh)

        # Prevent double output if the user has also configured the root logger
        logger.propagate = False

    return logger
