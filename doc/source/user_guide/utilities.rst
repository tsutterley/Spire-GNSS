============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Checks ``MD5`` or ``sha1`` hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/spire_toolkit/utilities.py


General Methods
===============

.. autofunction:: spire_toolkit.utilities.get_data_path

.. autofunction:: spire_toolkit.utilities.get_hash

.. autofunction:: spire_toolkit.utilities.url_split

.. autofunction:: spire_toolkit.utilities.get_unix_time

.. autofunction:: spire_toolkit.utilities.isoformat

.. autofunction:: spire_toolkit.utilities.even

.. autofunction:: spire_toolkit.utilities.ceil

.. autofunction:: spire_toolkit.utilities.copy

.. autofunction:: spire_toolkit.utilities.check_ftp_connection

.. autofunction:: spire_toolkit.utilities.ftp_list

.. autofunction:: spire_toolkit.utilities.from_ftp

.. autofunction:: spire_toolkit.utilities.check_connection

.. autofunction:: spire_toolkit.utilities.http_list

.. autofunction:: spire_toolkit.utilities.from_http

.. autofunction:: spire_toolkit.utilities.build_opener

.. autofunction:: spire_toolkit.utilities.gmao_list
