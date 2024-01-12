import helper


def convert_symbols(symbols: list, hgnc_dump):
    """ Convert list of gene symbols into HGNC ids

    Args:
        symbols (list): List of gene symbols to convert
        hgnc_dump (pd.Dataframe): HGNC dataframe

    Returns:
        list: List for which each element is composed of gene symbol, HGNC id
        and status
    """

    hgnc_ids = []

    for symbol in symbols:
        hgnc_id = helper.find_hgnc_id(symbol, hgnc_dump)
        hgnc_ids.append(hgnc_id)

    return hgnc_ids

def symbols_to_hgcn_ids(symbols: list, hgnc_dump_path: str):
    hgnc_dump = helper.parse_tsv(hgnc_dump_path)
    symbols = [symbol.upper() for symbol in symbols]
    hgnc_ids = convert_symbols(symbols, hgnc_dump)
    return hgnc_ids