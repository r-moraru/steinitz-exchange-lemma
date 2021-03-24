def echelon_form(matrix):
    lead = 0
    rows = len(matrix)
    cols = len(matrix[0])

    for r in range(rows):
        if lead >= cols:
            break

        i = r

        while matrix[i][lead] == 0:
            i += 1
            if i == rows:
                i = r
                lead += 1
                if lead >= cols:
                    break

        matrix[i], matrix[r] = matrix[r], matrix[i]
        div = matrix[r][lead]

        matrix[r] = [element / div for element in matrix[r]]

        for j in range(rows):
            if j != r:
                mul = matrix[j][lead]
                matrix[j] = [elj - elr * mul for elj, elr in zip(matrix[j], matrix[r])]
        
        lead += 1


def multiply(left_matrix, right_matrix):
    n = len(left_matrix)
    t = len(left_matrix[0])
    m = len(right_matrix[0])

    result = []

    for i in range(n):
        a = []
        for j in range(m):
            s = 0
            for k in range(t):
                s += left_matrix[i][k] * right_matrix[k][j]
            
            a.append(s)
        result.append(a)

    return result


def erase_line(matrix, line, non_zero):
    for col in range(len(matrix[0])):
        if matrix[line][col] != 0:
            non_zero[col] -= 1
        matrix[line][col] = 0


def erase_col(matrix, col):
    for line in range(len(matrix)):
        matrix[line][col] = 0


def eval_line(coefficient_matrix, col1, line, non_zero):
    for col2 in range(len(coefficient_matrix[0])):
        if col2 != col1 and coefficient_matrix[line][col2]:
            if non_zero[col2] <= 1:
                return False

    erase_line(coefficient_matrix, line, non_zero)
    erase_col(coefficient_matrix, col1)
    return True


def from_coefficient_matrix(coefficient_matrix, lis_size):
    """Using the coordinates of vectors in linearly independent set over
    the basis (remaining subset of generating set after first transformation),
    determines next vectors to be removed from generating set."""

    n = len(coefficient_matrix)
    m = len(coefficient_matrix[0])

    non_zero = []
    ordered_cols = []

    for col in range(m):
        s = 0
        for row in range(n):
            if coefficient_matrix[row][col] != 0:
                s += 1
        non_zero.append(s)
        ordered_cols.append([s, col])

    ordered_cols.sort(key=lambda x: x[0])

    exchangeable_lines = []

    for i in range(len(ordered_cols)):
        col1 = ordered_cols[i][1]
        for line in range(n):
            if coefficient_matrix[line][col1] != 0:
                if eval_line(coefficient_matrix, col1, line, non_zero):
                    exchangeable_lines.append(line)
                    break

    return exchangeable_lines


def get_transpose(matrix):
    rows = len(matrix)
    cols = len(matrix[0])

    transpose_matrix = [[matrix[j][i] for j in range(rows)] for i in range(cols)]
    
    return transpose_matrix


def no_pivot(matrix):
    """Find rows with no pivot in a matrix that is already
    in reduced row echelon form."""

    pivot_row = 0
    rows = len(matrix)
    cols = len(matrix[0])

    free_cols = []

    for i in range(cols):
        if pivot_row < rows and matrix[pivot_row][i] == 1:
            pivot_row += 1
        else:
            free_cols.append(i)

    return free_cols


def move_vectors(matrix, vector_indices, max_moves):
    """Returns a set containing vectors found at given indices in matrix.
    Deletes those vectors from matrix."""

    max_dif = min(len(vector_indices), max_moves)
    vector_indices.sort(reverse=True)
    A = [matrix.pop(vector_indices[i]) for i in range(max_dif)]
    return A


def inverse(matrix):
    lead = 0
    m = len(matrix)
    n = len(matrix[0])

    inverse_matrix = [[matrix[i][j] for j in range(n)] for i in range(m)]

    for _ in range(n):
        line = [1 if j == lead else 0 for j in range(n)]
        lead += 1
        inverse_matrix.append(line)

    temp_transpose = get_transpose(inverse_matrix)
    echelon_form(temp_transpose)

    inverse_matrix = [[temp_transpose[i][j] for j in range(m, 2 * m)] for i in range(n)]

    return inverse_matrix
