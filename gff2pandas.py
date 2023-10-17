"""
    This is a module taken from:
    https://github.com/foerstner-lab/gffpandas/releases/tag/v1.2.2
    Version: v1.2.2

    The module contains a class Gff3DataFrame, which is used to read a gff3 file.
    The main function used is converting the attributes column into
    separate columns for a pandas dataframe.

    Changes made:
    _split_atts function modified with try/except to handle AttributeError.
    Line 64: dtype of start and end changed to np.uint32
        dtype={
        "seq_id": str,
        "source": "category",
        "type": "category",
        "start": np.uint32,
        "end": np.uint32,
        "score": str,
        "strand": str,
        "phase": str,
        "attributes": str,
        }
"""

import itertools
import pandas as pd
import numpy as np


def read_gff3(input_file):
    return Gff3DataFrame(input_file)


def _split_atts(atts):
    """Split a feature string into attributes."""
    if not atts:
        return {}
    try:
        splits_list = [a.split("=") for a in atts.split(";") if "=" in a]
        return {item[0]: "=".join(item[1:]) for item in splits_list}
    except AttributeError:
        # Handle the case when atts has no 'split' attribute
        return {}


class Gff3DataFrame(object):
    """
    This class contains header information in the header attribute and
    a actual annotation data in the pandas dataframe in the df
    attribute.
    """

    def __init__(
        self, input_gff_file=None, input_df=None, input_header=None
    ) -> None:
        """Create an instance."""
        if input_gff_file is not None:
            self._gff_file = input_gff_file
            self._read_gff3_to_df()
            self._read_gff_header()
        else:
            self.df = input_df
            self.header = input_header

    def _read_gff3_to_df(self) -> pd.DataFrame:
        """Create a pandas dataframe.

        By the pandas library the gff3 file is read and
        a pd dataframe with the given column-names is returned."""
        self.df = pd.read_csv(
            self._gff_file,
            sep="\t",
            comment="#",
            names=[
                "seq_id",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
            dtype={
                "seq_id": str,
                "source": "category",
                "type": "category",
                "start": np.uint32,
                "end": np.uint32,
                "score": str,
                "strand": str,
                "phase": str,
                "attributes": str,
            },
        )
        return self.df

    def _read_gff_header(self) -> str:
        """
        Create a header string.

        The header of the gff file is read, means all lines,
        which start with '#'."""
        self.header = ""
        with open(self._gff_file) as file:
            for line in file:
                if line.startswith("#"):
                    self.header += line
                else:
                    break  # If no '#' is found, the header is finished.
        return self.header

    def attributes_to_columns(self) -> pd.DataFrame:
        """Saving each attribute-tag to a single column.

        Attribute column will be split by the tags in the single columns.
        For this method only a pandas DataFrame and not a Gff3DataFrame
        will be returned. Therefore, this data frame can not be saved as
        gff3 file.
        at_dic: dictionary with the attribute tags as keys
        and the corrosponding values.

        :return: pandas dataframe, whereby the attribute column of the gff3
                file are splitted into the different attribute tags
        :rtype: pandas DataFrame
        """
        attribute_df = self.df.copy()
        df_attributes = attribute_df.loc[:, "seq_id":"attributes"]
        attribute_df["at_dic"] = attribute_df.attributes.apply(_split_atts)
        attribute_df["at_dic_keys"] = attribute_df["at_dic"].apply(
            lambda at_dic: list(at_dic.keys())
        )
        merged_attribute_list = list(
            itertools.chain.from_iterable(attribute_df["at_dic_keys"])
        )
        nonredundant_list = sorted(list(set(merged_attribute_list)))
        for atr in nonredundant_list:
            df_attributes[atr] = attribute_df["at_dic"].apply(
                lambda at_dic: at_dic.get(atr)
            )
        return df_attributes

    def stats_dic(self) -> dict:
        """Gives the following statistics for the data:

        The maximal bp-length, minimal bp-length, the count of sense (+) and
        antisense (-) strands as well as the count of each available
        feature.

        :return: information about the given dataframe, which are the length of
            the longest and shortest feature entry (in bp), the number of
            feature on the sense and antisense strand and the number of
            different feature types.
        :rtype: dictionary
        """
        df_w_region = self.df[self.df.type != "region"]
        gene_length = df_w_region.end - df_w_region.start
        strand_counts = pd.value_counts(self.df["strand"]).to_dict()
        type_counts = pd.value_counts(self.df["type"]).to_dict()
        stats_dic = {
            "Maximal_bp_length": gene_length.max(),
            "Minimal_bp_length": gene_length.min(),
            "Counted_strands": strand_counts,
            "Counted_feature_types": type_counts,
        }
        return stats_dic
