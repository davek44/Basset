import sys


def message(msg, msgtype='progress'):
    '''Send a message to console.
    Args:
        msgtype (str): one of 'progress', 'warning', 'error', or 'debug'
    '''
    message = '[%(level)s]: %(text)s' % dict(level=msgtype.upper(), text=msg)
    sys.stderr.write(message.strip() + '\n')

    if msgtype == 'error':
        sys.exit(1)
