import unittest


class BaseTest(unittest.TestCase):
    def setUp(self) -> None:
        super(BaseTest, self).setUp()
        print(f"Running test: {self._testMethodName}\n")
        pass
