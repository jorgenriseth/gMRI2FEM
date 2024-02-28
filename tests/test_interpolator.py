import subprocess


def test_interpolator():
    subprocess.run(
        "OMP_NUM_THREADS=1 mpirun -n 2 python tests/test_resources/interpolator_script.py", shell=True
    ).check_returncode()
