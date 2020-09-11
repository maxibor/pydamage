from pypmml import Model
import pkg_resources


def load_model():
    """Returns the pmml model"""
    # This is a stream,like object. If you want the actual info, call
    # stream.read()
    stream = pkg_resources.resource_stream(__name__, "models/accuracy_model.xml")
    return Model.load(stream)


def prepare_data(pd_df):
    """Prepare pydamage result data for accuracy modelling

    Args:
        pd_df (pandas DataFrame):pydamage df result
    """
    reflen_bins = [
        (500, 1000),
        (1000, 2000),
        (2000, 5000),
        (5000, 10000),
        (10000, 20000),
        (20000, 50000),
        (50000, 100000),
        (100000, 200000),
        (200000, 500000),
    ]

    coverage_bins = [
        (1, 2),
        (2, 3),
        (3, 5),
        (5, 10),
        (10, 20),
        (20, 50),
        (50, 100),
        (100, 200),
        (200, 500),
    ]
    pd_df = pd_df[["coverage", "reflen", "pmax", "gc_content"]]
