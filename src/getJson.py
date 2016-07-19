from urllib.request import urlopen
from urllib.parse import urlencode
import httplib2
import json


def return_json(url='http://150.163.55.239:8081/custom/wisajson'):
    response = urlopen(url).read().decode('utf-8')

    jsonData = json.loads(response)
    return jsonData

def return_json_http():
    h = httplib2.Http()
    response, content = h.request('http://150.163.55.239:8081/custom/wisajson')
    return response, content

def get_json():
    pass

if __name__ == '__main__':
