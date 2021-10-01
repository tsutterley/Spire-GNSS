gmao_spire_gnss_sync.py
=======================

- Syncs Spire GNSS grazing angle altimetry data from the [NASA Global Modeling and Assimilation Office (GMAO)](https://gmao.gsfc.nasa.gov)

#### Calling Sequence
```bash
python gmao_spire_gnss_sync.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/Spire-GNSS/blob/main/scripts/gmao_spire_gnss_sync.py)

#### Command Line Options
- `-U X`, `--user X`: username for NASA GMAO Extranet Login
- `-P X`, `--password X`: Password for NASA GMAO Extranet Login
- `-N X`, `--netrc X`: path to .netrc file for authentication
- `-D X`, `--directory X`: working data directory
- `-t X`, `--timeout X`: Timeout in seconds for blocking operations
- `-l`, `--log`: output log of files downloaded
- `-M X`, `--mode X`: permissions mode of the directories and files synced
