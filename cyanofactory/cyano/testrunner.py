from django.conf import settings
from django.test.runner import DiscoverRunner

# http://tinystruggles.com/2014/10/20/django_tests_with_mirror_database.html

class MirrorTestDBTestSuiteRunner(DiscoverRunner):
    def run_tests(self, test_labels, extra_tests=None, **kwargs):
        # This test runner tried to run test on everything,
        # even on the django itself, this is a hacky workaround ;)
        #super(MirrorTestDBTestSuiteRunner, self).run_tests(
        #    settings.TESTED_APPS, extra_tests, **kwargs)
        self.setup_test_environment()
        suite = self.build_suite(test_labels, extra_tests)
        #old_config = self.setup_databases()
        result = self.run_suite(suite)
        #self.teardown_databases(old_config)
        self.teardown_test_environment()
        return self.suite_result(suite, result)