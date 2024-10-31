"""Some random utility functions.
"""
import os
import requests
import urllib.parse
import pandas as pd


def is_url(path):
    try:
        result = urllib.parse.urlparse(path)
        return all([result.scheme, result.netloc])
    except ValueError:
        return False


def get_file_size(file_url):
    """Calculates the file size in MB.

    Args:
        file_url (str): URL or local file path for which the file size is to be computed.

    Returns:
        float: Returns the file size in MB.
    """
    if is_url(file_url):
        # Get the file size from the URL
        response = requests.head(file_url)
        file_size = int(response.headers.get('content-length', 0))
    
    else:   # file_url is local filepath
        # Get the file size in bytes
        file_size = os.path.getsize(file_url)

    # Convert the file size to megabytes
    file_size_mb = file_size / (1024 * 1024)

    return file_size_mb
