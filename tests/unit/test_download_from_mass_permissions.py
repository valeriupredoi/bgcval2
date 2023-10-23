import os
import shutil
import stat

from bgcval2.download_from_mass import download_from_mass as dfm


def test_allowed_user():
    user_home = os.path.expanduser('~')
    output_folder = "bgcval2/local_test/BGC_data/u-xxx/"
    run_dir = os.path.join(user_home, output_folder)
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    dfm(jobID="u-xxx",
        doMoo=False,
        auto_download=False,
        config_user="defaults"
    )
    assert os.stat(run_dir).st_mode == 16893


def test_disallowed_user():
    user_home = os.path.expanduser('~')
    output_folder = "bgcval2/local_test/BGC_data/u-yyy/"
    run_dir = os.path.join(user_home, output_folder)
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    st = os.stat(run_dir)
    os.chmod(run_dir, False)
    dfm(jobID="u-yyy",
        doMoo=False,
        auto_download=False,
        config_user="defaults"
    )
    assert os.stat(run_dir).st_mode == 16384
    shutil.rmtree(run_dir, ignore_errors=True)
