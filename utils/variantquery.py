import requests


def restquery(url, data=None, method="GET"):
    """
    Send request to RESTful server and return the response from server
    :param url: query
    :param data: if method is "POST", data in JSON array type should be provided
    :param method: "GET" or "POST", default: "GET"
    :return: response from server with JSON type
    """
    if method == "GET":
        req = requests.get(url)
    else:
        req = requests.post(url, data=data)
    return req.json()


def createquery(baseurl, ch, pos, ref, alt):
    """
    Create query with baseurl, chromosome number, position, reference allele, alternate allele
    :param baseurl:
    :param ch:
    :param pos:
    :param ref:
    :param alt:
    :return:
    """
    args = [ch, pos, ref, alt]
    query = '-'.join(args)
    url = '/'.join([baseurl, query])
    return url
