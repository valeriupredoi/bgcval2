import subprocess


def bgcval2_test_data():
    """reate test data (the big one)."""
    bash_command = "python tests/integration/produce_dummy_data.py"
    print(f"Creating dummy data by running the command: {bash_command}")
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output = str(process.communicate()[0])
    print(f"{output}")

    return output
