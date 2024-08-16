from math import sqrt
from sklearn.metrics import mean_squared_error


def dist(p1, p2):
    """2D Euclidean distance function

    Args:
        p1 (tuple): first point
        p2 (tuple): second point

    Returns:
        float: Euclidean distance
    """
    return sqrt(((p1[0]-p2[0])**2)+((p1[1]-p2[1])**2))


def mse(y_true, y_pred):
    """Error Metric: Mean Squared Error 

    Args:
        y_true (numpy.array): true values
        y_pred (numpy.array): predicted values

    Returns:
        float: mean squared error
    """
    return mean_squared_error(y_true, y_pred)