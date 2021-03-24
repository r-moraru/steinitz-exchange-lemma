import helpers


def find_replaceable_subset(gs, lis):
    """Given a generating set (gs) and a linearly independent set of vectors (lis)
    returns subset of generating set (exchangeable subset) such that
    lis U (gs - exchangeable_subset) generates the same vector space and has
    the same cardinality as gs"""

    lis_size = len(lis)
    exchangeable_subset = []

    # gs_matrix is a matrix containing vectors of the generating set placed vertically
    gs_matrix = helpers.get_transpose(gs)
    helpers.echelon_form(gs_matrix)

    # Find columns with no pivot
    free_cols = helpers.no_pivot(gs_matrix)

    # A base is obtained by removing all redundant vectors from a generating set
    # TODO: improve explanation
    exchangeable_subset += helpers.move_vectors(gs, free_cols, lis_size)

    # Number of vectors to be removed from linearly independent set
    num_remove = min(len(free_cols), lis_size)

    for i in range(num_remove-1, -1, -1):
        lis.pop(i)

    if len(lis) != 0:
        # Finding solution for linear equation: gs*coefficient_matrix = lis
        inv_gs_matrix = helpers.inverse(gs)
        coefficient_matrix = helpers.multiply(inv_gs_matrix, helpers.get_transpose(lis))

        # Find vectors to be exchanged
        remaining_vectors = helpers.from_coefficient_matrix(coefficient_matrix, len(lis))

        # Move remaining exchangeable vectors from generating set to exchangeable_subset
        exchangeable_subset += helpers.move_vectors(gs, remaining_vectors, len(lis))
    
    return exchangeable_subset
