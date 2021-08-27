import subprocess

TESTS = {
    '1 jacobi tests/input_J_0.txt': 'tests/output_J_0.txt',
    '1 jacobi tests/input_J_1.txt': 'tests/output_J_1.txt',
    '1 jacobi tests/input_J_2.txt': 'tests/output_J_2.txt',
    '1 jacobi tests/input_J_3.txt': 'tests/output_J_3.txt',
    '1 jacobi tests/input_J_4.txt': 'tests/output_J_4.txt',
    '1 jacobi tests/input_J_5.txt': 'tests/output_J_5.txt',
    '1 jacobi tests/input_J_6.txt': 'tests/output_J_6.txt',
    '1 jacobi tests/input_J_7.txt': 'tests/output_J_7.txt',
    '1 jacobi tests/input_J_8.txt': 'tests/output_J_8.txt',
    '1 jacobi tests/input_J_9.txt': 'tests/output_J_9.txt',
}


def _run_process(cmd):
    print('> Running \'{}\'... '.format(cmd), end='')
    proc = subprocess.run(cmd, shell=True, capture_output=True)
    print('done')

    if proc.returncode != 0:
        raise RuntimeError('Invalid return code for command line \'{}\': {}'.format(
            cmd, proc.returncode))
    return proc.stdout.decode('utf8')


def _read_file(path):
    with open(path, 'r') as fd:
        return fd.read()


def main():
    success = True
    for exec_command in ['python3 spkmeans.py', './spkmeans']:
        for args, expected_output_file in TESTS.items():
            full_command = '{} {}'.format(exec_command, args)
            output = _run_process(full_command)
            expected_output = _read_file(expected_output_file)
            if output != expected_output:
                success = False
                print('\nMismatching outputs for command line \'{}\':'.format(full_command))
                print('Expected:\n{}'.format(expected_output))
                print('Got:\n{}'.format(output))
    print('>>> {} <<<'.format('Success' if success else 'Failure'))


if __name__ == '__main__':
    main()
