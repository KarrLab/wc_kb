""" Configuration

:Author: Jonathan Karr <onrkarr@gmail.com>
:Date: 2019-01-06
:Copyright: 2019, Karr Lab
:License: MIT
"""

import configobj
import os
import pkg_resources
import wc_utils.config


def get_config(extra=None):
    """ Get configuration

    Args:
        extra (:obj:`dict`, optional): additional configuration to override

    Returns:
        :obj:`configobj.ConfigObj`: nested dictionary with the configuration settings loaded from the configuration source(s).
    """
    paths = wc_utils.config.ConfigPaths(
        default=pkg_resources.resource_filename('wc_kb', 'config/core.default.cfg'),
        schema=pkg_resources.resource_filename('wc_kb', 'config/core.schema.cfg'),
        user=(
            'wc_kb.cfg',
            os.path.expanduser('~/.wc/wc_kb.cfg'),
        ),
    )

    config = wc_utils.config.ConfigManager(paths).get_config(extra=extra)
    return config
