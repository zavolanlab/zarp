# This script targets the client api version 0.4.0 and later

#
#  Check the page: https://github.com/LabKey/labkey-api-python/blob/master/samples/query_examples.py
#  for example about filtering in queries.
#  A starting point to investigate further is here:
#  https://www.labkey.org/download/clientapi_docs/javascript-api/symbols/LABKEY.Query.Filter.html

import labkey
import pandas as pd
import sys

# for convenience, load QueryFilter explicitly (avoids long lines in filter definitions)
from labkey.query import QueryFilter

if __name__ == "__main__":
  # These are values of variables for which the script works
  # project_name = "TEST_ABOERSCH"
  # query_name = "RNA_Seq_data_template"
  project_name = sys.argv[1]
  query_name = sys.argv[2]
  server_context = labkey.utils.create_server_context('labkey.scicore.unibas.ch', '/Zavolan Group/'+project_name, 'labkey', use_ssl=True)
  schema_name = "lists"
  results = labkey.query.select_rows(server_context,schema_name,query_name)
  table_of_data = pd.DataFrame(results["rows"])
  print(table_of_data)

