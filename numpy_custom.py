import numpy as np


def gradient(f, res=1, outside='none'):
    # Wenn Eingabe fehlt: Outside Gebiete wo Array = 0
    if outside is 'none':
        outside = f == 0

    # Array Size
    fshape = f.shape
    # Gradient normal bestimmen
    grad_f = np.gradient(f, res)
    grad_f_0 = grad_f[0]
    grad_f_1 = grad_f[1]
    # Gradient ausserhalb auf 0 setzen
    grad_f_0[outside] = 0
    grad_f_1[outside] = 0
    # Differenz entlang Achse 0 bestimmen und durch res teilen, erweitern auf Arraygroesse, links - und rechtsseitig
    diff_f_0 = np.diff(f, axis=0) / res
    diff_f_0_left = np.append(diff_f_0, np.zeros((1, f.shape[1])), axis=0)
    diff_f_0_right = np.append(np.zeros((1, f.shape[1])), diff_f_0, axis=0)
    # differenz entlang Achse 1 bestimmen
    diff_f_1 = np.diff(f, axis=1) / res
    diff_f_1_left = np.append(diff_f_1, np.zeros((f.shape[0], 1)), axis=1)
    diff_f_1_right = np.append(np.zeros((f.shape[0], 1)), diff_f_1, axis=1)

    # Mithilfe diff den Rand des Gebiets ermitteln, fuer links-/rechtsseitig erweitern
    diff_outside_0 = np.diff(outside, axis=0)
    diff_outside_0_left = np.append(np.ones((1, f.shape[1]), dtype=bool), diff_outside_0, axis=0)
    diff_outside_0_right = np.append(diff_outside_0, np.ones((1, fshape[1]), dtype=bool), axis=0)
    diff_outside_1 = np.diff(outside, axis=1)
    diff_outside_1_left = np.append(np.ones((f.shape[0], 1), dtype=bool), diff_outside_1, axis=1)
    diff_outside_1_right = np.append(diff_outside_1, np.ones((fshape[0], 1), dtype=bool), axis=1)

    # Indizes links/rechtsseitige Raender bestimmen
    left_replace_0 = np.logical_and(diff_outside_0_left, np.logical_not(outside))
    right_replace_0 = np.logical_and(diff_outside_0_right, np.logical_not(outside))

    left_replace_1 = np.logical_and(diff_outside_1_left, np.logical_not(outside))
    right_replace_1 = np.logical_and(diff_outside_1_right, np.logical_not(outside))

    # Standardgradient an Raendern durch links/rechtsseitige Quotienten ersetzen
    grad_f_0[left_replace_0] = diff_f_0_left[left_replace_0]
    grad_f_0[right_replace_0] = diff_f_0_right[right_replace_0]

    grad_f_1[left_replace_1] = diff_f_1_left[left_replace_1]
    grad_f_1[right_replace_1] = diff_f_1_right[right_replace_1]

    # Output analog zu Numpy Funktion np.gradient
    return grad_f_0, grad_f_1


def cholesky(r_ij):
    a_ij = np.zeros(r_ij.shape)

    for i in range(a_ij.shape[0]):
        for j in range(i + 1):
            summe = r_ij[i, j]
            for k in range(j):
                summe = summe - a_ij[i, k] * a_ij[j, k]
            if i > j and a_ij[j, j] > 0:
                a_ij[i, j] = summe / a_ij[j, j]
            elif summe > 0 and i == j:
                a_ij[i, i] = np.sqrt(summe)

    return a_ij


# Aus Wikipedia Pseudocode fuer Cholesky Zerlegung
#     For i = 1 To n
#         For j = 1 To i
#             Summe = a(i, j)
#             For k = 1 To j-1
#                 Summe = Summe - a(i, k) * a(j, k)
#             If i > j Then
#                 a(i, j) = Summe / a(j, j)   // Untere Dreiecksmatrix
#             Else If Summe > 0 Then          // Diagonalelement
#                 a(i, i) = Sqrt(Summe)       // ... ist immer groesser Null
#             Else
#                 ERROR


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))