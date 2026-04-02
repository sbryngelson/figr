#!/usr/bin/env python3

import signal

# Only import what's needed for startup - other modules are loaded lazily
from figr import args, lock, state
from figr.common import FigrException, does_command_exist, quit, setup_debug_logging
from figr.printer import cons
from figr.state import ARG


def __print_greeting():
    pass


def __checks():
    if not does_command_exist("cmake"):
        raise FigrException("CMake is required to build MFC but couldn't be located on your system. Please ensure it installed and discoverable (e.g in your system's $PATH).")


def __run():
    # Lazy import modules only when needed for the specific command
    cmd = ARG("command")

    if cmd == "test":
        from figr.test import test

        test.test()
    elif cmd == "run":
        from figr.run import run

        run.run()
    elif cmd == "build":
        from figr import build

        build.build()
    elif cmd == "clean":
        from figr import clean

        clean.clean()
    elif cmd == "validate":
        from figr import validate

        validate.validate()


if __name__ == "__main__":
    try:
        lock.init()
        state.gARG = args.parse(state.gCFG)

        # Setup debug logging if requested
        setup_debug_logging(ARG("debug_log", dflt=None))

        # --reldebug and --debug are mutually exclusive: if the user explicitly
        # passed one, clear the other (which may be lingering from a prior run's
        # persisted config).
        if state.gARG.get("reldebug") and state.gARG.get("debug"):
            # Determine which flag the user explicitly passed on this invocation
            # by checking whether the persisted config already had it set.
            if not state.gCFG.reldebug:
                # User just passed --reldebug; clear persisted --debug
                state.gARG["debug"] = False
            else:
                # User just passed --debug; clear persisted --reldebug
                state.gARG["reldebug"] = False

        lock.switch(state.MFCConfig.from_dict(state.gARG))

        pass

        __print_greeting()
        __checks()
        __run()

    except FigrException as exc:
        cons.reset()
        cons.print(f"""\


[bold red]Error[/bold red]: {str(exc)}
""")
        quit(signal.SIGTERM)
    except KeyboardInterrupt:
        quit(signal.SIGTERM)
    except Exception as exc:
        cons.reset()
        cons.print_exception()
        cons.print(f"""\


[bold red]ERROR[/bold red]: An unexpected exception occurred: {str(exc)}
""")

        quit(signal.SIGTERM)
