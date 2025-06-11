def wmfd(tree1, tree2, weights):
    """
    Compute the Weighted Multi-Feature Distance (WMFD) between two lineage trees.

    Parameters:
        tree1 (dict): feature values of tree 1 (e.g., {'topology': 0.5, 'length': 0.3, ...})
        tree2 (dict): feature values of tree 2
        weights (dict): weights for each feature (e.g., {'topology': 0.6, 'length': 0.4, ...})

    Returns:
        float: the weighted sum of feature-wise differences
    """
    total = 0.0
    for feature in weights:
        if feature in tree1 and feature in tree2:
            diff = abs(tree1[feature] - tree2[feature])
            total += weights[feature] * diff
        else:
            raise ValueError(f"Missing feature '{feature}' in one of the trees.")
    return total
