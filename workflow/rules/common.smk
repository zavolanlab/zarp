## Function definitions
def get_sample(column_id, search_id=None, search_value=None):
    """Get relevant per sample information from samples table"""
    if search_id:
        if search_id == "index":
            return str(
                samples_table.loc[samples_table.index == search_value, column_id].iloc[
                    0
                ]
            )
        else:
            return str(
                samples_table.loc[
                    samples_table[search_id] == search_value, column_id
                ].iloc[0]
            )
    else:
        return str(samples_table.loc[0, column_id])


def get_directionality(libtype, tool):
    """Get directionality value for different tools"""
    directionality = ""

    for key in directionality_dict.keys():
        # Use the first of 'SF' or 'SR' that is found in libtype to look up
        # directionality params for the current tool
        if key in libtype:
            directionality = directionality_dict[key][tool]
            break

    # If libtype contains neither 'SF', nor 'SR' we don't know what to do
    if not directionality:
        raise ValueError(
            f"Unknown libtype {libtype}.\n"
            "Provide one of 'SF', 'SR', 'ISF', 'ISR', 'OSF', 'OSR', 'MSF', 'MSR' in samples.tsv!"
        )

    return directionality
