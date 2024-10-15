import argparse
import getpass
from COMPS import Client
from COMPS.CredentialPrompt import CredentialPrompt
import simulations.manifest as manifest
__comps_client_version = 10


class StaticCredentialPrompt(CredentialPrompt):
    def __init__(self, comps_url, username, password):
        if (comps_url is None) or (username is None) or (password is None):
            print("Usage: python create_auth_token_args.py --comps_url url --username username --password pwd")
            print("\n")
            raise RuntimeError('Missing comps_url, or username or password')
        self._times_prompted = 0
        self.comps_url = comps_url
        self.username = username
        self.password = password

    def prompt(self):
        print("logging in with hardcoded user/pw")
        self._times_prompted = self._times_prompted + 1
        if self._times_prompted > 3:
            raise RuntimeError('Failure authenticating')
        print("Hit here")
        return {'Username': self.username, 'Password': self.password}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--comps_url', '-c', default=manifest.endpoint, help='comps url')
    parser.add_argument('--username', '-u', help='enter username')
    parser.add_argument('--password', '-p', help='enter password')

    args = parser.parse_args()
    if args.password is None:
        p = getpass.getpass()
    else:
        p = args.password

    compshost = args.comps_url
    print(f"Login to {compshost} and cache the credentials.")

    Client.login(compshost, StaticCredentialPrompt(comps_url=args.comps_url, username=args.username,
                                                   password=p))
