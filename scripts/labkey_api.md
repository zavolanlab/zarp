> **NOTE**: Include this info in the main doc in the root directory once 
> available; add information on how to get credentials for the LabKey server 
> (i.e., how to obatain a password for the `.netrc` file)

In order to connect to the LabKey through API, you will first need to create a
file `.netrc` in your home directory:

```bash
touch ${HOME}/.netrc
```

Add the following lines to the file:

```console
machine <remote-instance-of-labkey-server>  
login <user-email>
password <user-password>  
```

To secure the file, set permissions in a way that only you can see the content
of the file: 

```bash
chmod 400 .netrc
```

Install the `labkey` and `pandas` packages, ideally from a virtual environment
(e.g., `virtualenv` or `conda`):

```bash
pip install labkey pandas
```

Run the LabKey API client script:

```bash
python labkey_api.py project_name labkey_table_nane  
```

Example:

```bash
python labkey_api.py TEST_ABOERSCH RNA_Seq_data_template  
```

Right now the script prints a representation of a `pandas` data frame
containing the requested LabKey table the the screen. For further processing
the current script could be included in another script, or it could be modified
to write out the data in a desired file format (e.g., TSV).
