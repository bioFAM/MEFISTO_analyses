# /opt/anaconda3/envs/qiime2-2020.8/bin

# The preprocessing script follows the preprocessing conducted by Martino et al (Nat Biotech. 2020) on the ECAM data.
# The input data and original preprocessing scripts are stored at ‘Code Ocean’ capsule: https://doi.org/10.24433/CO.5938114.v1.
# by Martino et al.

import biom
import skbio
from pandas import concat
from pandas import DataFrame
from skbio import OrdinationResults, DistanceMatrix
import qiime2 as q2
from gemelli.factorization import TensorFactorization
from gemelli.preprocessing import build, rclr
from gemelli._ctf_defaults import (DEFAULT_COMP, DEFAULT_MSC,
                                   DEFAULT_MFC, DEFAULT_MAXITER,
                                   DEFAULT_FMETA as DEFFM)
import pickle


def preprocess(table: biom.Table,
        sample_metadata: DataFrame,
        individual_id_column: str,
        state_column: str,
        min_sample_count: int = DEFAULT_MSC,
        min_feature_count: int = DEFAULT_MFC,
        feature_metadata: DataFrame = DEFFM):

    state_columns = [state_column]

    # transform artifact to biom table
    table = table.view(biom.Table)


    # validate the metadata using q2 as a wrapper
    if sample_metadata is not None and not isinstance(sample_metadata,
                                                      DataFrame):
        sample_metadata = sample_metadata.to_dataframe()
    keep_cols = state_columns + [individual_id_column]
    all_sample_metadata = sample_metadata.drop(keep_cols, axis=1)
    sample_metadata = sample_metadata[keep_cols]
    # validate the metadata using q2 as a wrapper
    if feature_metadata is not None and not isinstance(feature_metadata,
                                                       DataFrame):
        feature_metadata = feature_metadata.to_dataframe()
    # match the data (borrowed in part from gneiss.util.match)
    subtablefids = table.ids('observation')
    subtablesids = table.ids('sample')
    if len(subtablesids) != len(set(subtablesids)):
        raise ValueError('Data-table contains duplicate sample IDs')
    if len(subtablefids) != len(set(subtablefids)):
        raise ValueError('Data-table contains duplicate feature IDs')
    submetadataids = set(sample_metadata.index)
    subtablesids = set(subtablesids)
    subtablefids = set(subtablefids)
    if feature_metadata is not None:
        submetadatafeat = set(feature_metadata.index)
        fidx = subtablefids & submetadatafeat
        if len(fidx) == 0:
            raise ValueError(("No more features left.  Check to make "
                              "sure that the sample names between "
                              "`feature-metadata` and `table` are "
                              "consistent"))
        feature_metadata = feature_metadata.reindex(fidx)
    sidx = subtablesids & submetadataids
    if len(sidx) == 0:
        raise ValueError(("No more features left.  Check to make sure that "
                          "the sample names between `sample-metadata` and"
                          " `table` are consistent"))
    if feature_metadata is not None:
        table.filter(list(fidx), axis='observation', inplace=True)
    table.filter(list(sidx), axis='sample', inplace=True)
    sample_metadata = sample_metadata.reindex(sidx)

    # filter and import table
    for axis, min_sum in zip(['sample',
                              'observation'],
                             [min_sample_count,
                              min_feature_count]):
        table = table.filter(table.ids(axis)[table.sum(axis) >= min_sum],
                             axis=axis, inplace=True)

    # table to dataframe
    table = DataFrame(table.matrix_data.toarray(),
                      table.ids('observation'),
                      table.ids('sample'))

    pickle.dump(table, open( "data/processed_data/ECMA_table.p", "wb" ) )
    # tensor building
    tensor = build()
    tensor.construct(table, sample_metadata,
                     individual_id_column, state_columns)

    return(tensor)


# preprocess ECAM data
if __name__ == "__main__":
    datadir = "data/CTF-data/ECAM/"
    # import table
    table = q2.Artifact.load(datadir + 'table.qza')

    # import sample metadata
    metadata = q2.Metadata.load(datadir + 'metadata-matched.tsv')

    # import feature metadata
    taxonomy = q2.Artifact.load(datadir + 'taxonomy.qza')

    # save meta data
    sample_meta = metadata.to_dataframe()
    sample_meta.to_csv("data/processed_data/sample_metadata.csv" )
    feature_meta = taxonomy.view(q2.Metadata).to_dataframe()
    feature_meta.to_csv("data/processed_data/taxonomy.csv" )

    # build tensor
    tensor = preprocess(table,
               metadata,
               "host_subject_id",
               "month",
               min_sample_count =0,
               min_feature_count=5,
               feature_metadata=taxonomy.view(q2.Metadata)
               )
    pickle.dump(tensor, open( "data/processed_data/ECMA_tensor.p", "wb" ) )

    # do rclr transformation and save this for use in MEFISTO
    rclr_data = rclr(tensor.counts)
    print(rclr_data.shape)
    pickle.dump(rclr_data, open( "data/processed_data/ECAM_rclr_counts.p", "wb" ))
