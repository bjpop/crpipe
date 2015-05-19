'''
Retrieve the version number of the pipeline from the package
configuration.
'''

# This solution was suggested on Stack Overflow:
# http://stackoverflow.com/questions/2058802/how-can-i-get-the-version-defined-in-setup-py-setuptools-in-my-package

import pkg_resources  # part of setuptools

version = pkg_resources.require("crpipe")[0].version
