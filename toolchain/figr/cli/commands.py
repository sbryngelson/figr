"""
figr CLI Command Definitions - SINGLE SOURCE OF TRUTH

All command definitions live here. This file is used to generate:
- argparse parsers
- bash/zsh completions
- user guide help content
- CLI reference documentation

When adding a new command or option, ONLY modify this file.
Then run `./figr.sh generate` to update completions.
"""

from .schema import ArgAction, Argument, CLISchema, Command, CommonArgumentSet, Completion, CompletionType, Example, MutuallyExclusiveGroup, Positional

# CONSTANTS (shared with other modules)

TARGET_NAMES = ["hdf5", "silo", "pre_process", "simulation", "post_process", "syscheck", ]

DEFAULT_TARGET_NAMES = ["pre_process", "simulation", "post_process"]

TEMPLATE_NAMES = ["bridges2", "carpenter", "carpenter-cray", "default", "delta", "deltaai", "frontier", "hipergator", "nautilus", "oscar", "phoenix", "phoenix-bench", "santis", "tuo"]

GPU_OPTIONS = ["acc", "mp"]

ENGINE_OPTIONS = ["interactive", "batch"]

MPI_BINARIES = ["mpirun", "jsrun", "srun", "mpiexec"]


# COMMON ARGUMENT SETS

COMMON_TARGETS = CommonArgumentSet(
    name="targets",
    arguments=[
        Argument(
            name="targets",
            short="t",
            help="Space-separated list of targets to act upon.",
            nargs="+",
            type=str,
            default=DEFAULT_TARGET_NAMES,
            choices=TARGET_NAMES,
            metavar="TARGET",
            completion=Completion(type=CompletionType.CHOICES, choices=TARGET_NAMES),
        ),
    ],
)

COMMON_JOBS = CommonArgumentSet(
    name="jobs",
    arguments=[
        Argument(
            name="jobs",
            short="j",
            help="Allows for JOBS concurrent jobs.",
            type=int,
            default=1,
            metavar="JOBS",
        ),
    ],
)

COMMON_VERBOSE = CommonArgumentSet(
    name="verbose",
    arguments=[
        Argument(
            name="verbose",
            short="v",
            help="Increase output verbosity (-v, -vv, -vvv for more detail).",
            action=ArgAction.COUNT,
            default=0,
        ),
    ],
)

COMMON_DEBUG_LOG = CommonArgumentSet(
    name="debug_log",
    arguments=[
        Argument(
            name="debug-log",
            short="d",
            help="Enable Python toolchain debug logging (not figr code).",
            action=ArgAction.STORE_TRUE,
            dest="debug_log",
        ),
    ],
)

COMMON_GPUS = CommonArgumentSet(
    name="gpus",
    arguments=[
        Argument(
            name="gpus",
            short="g",
            help="(Optional GPU override) List of GPU #s to use (environment default if unspecified).",
            nargs="+",
            type=int,
            default=None,
        ),
    ],
)

# FigrConfig flags are handled specially in argparse_gen.py
# This marker tells the generator to add --mpi/--no-mpi, --gpu/--no-gpu, etc.
COMMON_FIGR_CONFIG = CommonArgumentSet(
    name="figr_config",
    figr_config_flags=True,
    arguments=[],  # Generated dynamically
)


# COMMAND DEFINITIONS

BUILD_COMMAND = Command(
    name="build",
    aliases=["b"],
    help="Build figr and its dependencies.",
    description="Build figr targets with optional GPU support.",
    include_common=["targets", "figr_config", "jobs", "verbose", "debug_log"],
    arguments=[
        Argument(
            name="deps-only",
            help="Only fetch and build dependencies, do not build figr targets.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="deps_only",
        ),
    ],
    examples=[
        Example("./figr.sh build", "Build all default targets (CPU)"),
        Example("./figr.sh build -j 8", "Build with 8 parallel jobs"),
        Example("./figr.sh build --gpu", "Build with GPU (OpenACC) support"),
    ],
    key_options=[
        ("-j, --jobs N", "Number of parallel build jobs"),
        ("-t, --targets", "Targets: pre_process, simulation, post_process"),
        ("--gpu [acc|mp]", "Enable GPU support (OpenACC or OpenMP)"),
        ("--debug", "Build in debug mode"),
    ],
)

RUN_COMMAND = Command(
    name="run",
    aliases=["r"],
    help="Run a case with figr.",
    description="Run a figr simulation case interactively or submit as a batch job.",
    include_common=["targets", "figr_config", "jobs", "verbose", "debug_log", "gpus"],
    positionals=[
        Positional(
            name="input",
            help="Input file to run.",
            completion=Completion(type=CompletionType.FILES_PY),
        ),
    ],
    arguments=[
        Argument(
            name="engine",
            short="e",
            help="Job execution/submission engine choice.",
            choices=ENGINE_OPTIONS,
            default="interactive",
            completion=Completion(type=CompletionType.CHOICES, choices=ENGINE_OPTIONS),
        ),
        Argument(
            name="partition",
            short="p",
            help="(Batch) Partition for job submission.",
            default="",
            metavar="PARTITION",
        ),
        Argument(
            name="quality_of_service",
            short="q",
            help="(Batch) Quality of Service for job submission.",
            default="",
            metavar="QOS",
        ),
        Argument(
            name="nodes",
            short="N",
            help="(Batch) Number of nodes.",
            type=int,
            default=1,
            metavar="NODES",
        ),
        Argument(
            name="tasks-per-node",
            short="n",
            help="Number of tasks per node.",
            type=int,
            default=1,
            metavar="TASKS",
            dest="tasks_per_node",
        ),
        Argument(
            name="walltime",
            short="w",
            help="(Batch) Walltime.",
            default="01:00:00",
            metavar="WALLTIME",
        ),
        Argument(
            name="account",
            short="a",
            help="(Batch) Account to charge.",
            default="",
            metavar="ACCOUNT",
        ),
        Argument(
            name="email",
            short="@",
            help="(Batch) Email for job notification.",
            default="",
            metavar="EMAIL",
        ),
        Argument(
            name="name",
            short="#",
            help="(Batch) Job name.",
            default="figr",
            metavar="NAME",
        ),
        Argument(
            name="scratch",
            short="s",
            help="Build from scratch.",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
        Argument(
            name="binary",
            short="b",
            help="(Interactive) Override MPI execution binary",
            choices=MPI_BINARIES,
            default=None,
            completion=Completion(type=CompletionType.CHOICES, choices=MPI_BINARIES),
        ),
        Argument(
            name="dry-run",
            help="(Batch) Run without submitting batch file.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="dry_run",
        ),
        Argument(
            name="no-build",
            help="(Testing) Do not rebuild figr.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="no_build",
        ),
        Argument(
            name="wait",
            help="(Batch) Wait for the job to finish.",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
        Argument(
            name="computer",
            short="c",
            help="(Batch) Path to a custom submission file template or one of the built-in templates.",
            default="default",
            metavar="COMPUTER",
            completion=Completion(type=CompletionType.CHOICES, choices=TEMPLATE_NAMES),
        ),
        Argument(
            name="output-summary",
            short="o",
            help="Output file (YAML) for summary.",
            default=None,
            metavar="OUTPUT",
            dest="output_summary",
        ),
        Argument(
            name="clean",
            help="Clean the case before running.",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
        # Profiler arguments with REMAINDER
        Argument(
            name="ncu",
            help="Profile with NVIDIA Nsight Compute.",
            nargs="...",  # REMAINDER
            type=str,
        ),
        Argument(
            name="nsys",
            help="Profile with NVIDIA Nsight Systems.",
            nargs="...",  # REMAINDER
            type=str,
        ),
        Argument(
            name="rcu",
            help="Profile with ROCM rocprof-compute.",
            nargs="...",  # REMAINDER
            type=str,
        ),
        Argument(
            name="rsys",
            help="Profile with ROCM rocprof-systems.",
            nargs="...",  # REMAINDER
            type=str,
        ),
    ],
    examples=[
        Example("./figr.sh run case.py", "Run interactively with 1 rank"),
        Example("./figr.sh run case.py -n 4", "Run with 4 MPI ranks"),
        Example("./figr.sh run case.py -e batch -N 2 -n 4", "Submit batch job: 2 nodes, 4 ranks/node"),
    ],
    key_options=[
        ("-n, --tasks-per-node", "MPI ranks per node"),
        ("-N, --nodes", "Number of nodes (batch)"),
        ("-e, --engine", "interactive or batch"),
        ("-a, --account", "Account to charge (batch)"),
        ("-w, --walltime", "Wall time limit (batch)"),
    ],
)

TEST_COMMAND = Command(
    name="test",
    aliases=["t"],
    help="Run figr's test suite.",
    description="Run figr's test suite with various filtering and generation options.",
    include_common=["figr_config", "jobs", "verbose", "debug_log", "gpus"],
    # Note: does NOT include "targets" - test uses different target handling
    arguments=[
        Argument(
            name="list",
            short="l",
            help="List all available tests.",
            action=ArgAction.STORE_TRUE,
        ),
        Argument(
            name="from",
            short="f",
            help="First test UUID to run.",
            default=None,
            type=str,
        ),
        Argument(
            name="to",
            short="t",
            help="Last test UUID to run.",
            default=None,
            type=str,
        ),
        Argument(
            name="only",
            short="o",
            help="Only run tests with specified properties.",
            nargs="+",
            type=str,
            default=[],
            metavar="L",
        ),
        Argument(
            name="test-all",
            short="a",
            help="Run the Post Process Tests too.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="test_all",
        ),
        Argument(
            name="percent",
            short="%",
            help="Percentage of tests to run.",
            type=int,
            default=100,
        ),
        Argument(
            name="max-attempts",
            short="m",
            help="Maximum number of attempts to run a test.",
            type=int,
            default=1,
            dest="max_attempts",
        ),
        Argument(
            name="rdma-mpi",
            help="Run tests with RDMA MPI enabled",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="rdma_mpi",
        ),
        Argument(
            name="no-build",
            help="(Testing) Do not rebuild figr.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="no_build",
        ),
        Argument(
            name="no-examples",
            help="Do not test example cases.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="no_examples",
        ),
        Argument(
            name="dry-run",
            help="Build and generate case files but do not run tests.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="dry_run",
        ),
        Argument(
            name="shard",
            help="Run only a subset of tests (e.g., '1/2' for first half, '2/2' for second half).",
            type=str,
            default=None,
        ),
        Argument(
            name="build-coverage-cache",
            help="Run all tests with gcov instrumentation to build the file-level coverage cache. Pass --gcov to enable coverage instrumentation in the internal build step.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="build_coverage_cache",
        ),
        Argument(
            name="only-changes",
            help="Only run tests whose covered files overlap with files changed since branching from master (uses file-level gcov coverage cache).",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="only_changes",
        ),
        Argument(
            name="changes-branch",
            help="Branch to compare against for --only-changes (default: master).",
            type=str,
            default="master",
            dest="changes_branch",
        ),
    ],
    mutually_exclusive=[
        MutuallyExclusiveGroup(
            arguments=[
                Argument(
                    name="generate",
                    help="(Test Generation) Generate golden files.",
                    action=ArgAction.STORE_TRUE,
                    default=False,
                ),
                Argument(
                    name="add-new-variables",
                    help="(Test Generation) If new variables are found in D/ when running tests, add them to the golden files.",
                    action=ArgAction.STORE_TRUE,
                    default=False,
                    dest="add_new_variables",
                ),
                Argument(
                    name="remove-old-tests",
                    help="(Test Generation) Delete test directories that are no longer needed.",
                    action=ArgAction.STORE_TRUE,
                    default=False,
                    dest="remove_old_tests",
                ),
            ]
        ),
    ],
    examples=[
        Example("./figr.sh test", "Run all tests"),
        Example("./figr.sh test -j 4", "Run with 4 parallel jobs"),
        Example("./figr.sh test --only 3D", "Run only 3D tests"),
        Example("./figr.sh test --generate", "Regenerate golden files"),
        Example("./figr.sh test --only-changes -j 4", "Run tests affected by changed files"),
        Example("./figr.sh build --gcov -j 8 && ./figr.sh test --build-coverage-cache", "One-time: build file-coverage cache"),
    ],
    key_options=[
        ("-j, --jobs N", "Number of parallel test jobs"),
        ("-o, --only PROP", "Run tests matching property"),
        ("-f, --from UUID", "Start from specific test"),
        ("--generate", "Generate/update golden files"),
        ("--no-build", "Skip rebuilding figr"),
        ("--build-coverage-cache", "Build file-level gcov coverage cache (one-time)"),
        ("--only-changes", "Run tests affected by changed files (requires cache)"),
    ],
)

CLEAN_COMMAND = Command(
    name="clean",
    aliases=["c"],
    help="Clean build artifacts.",
    description="Remove build artifacts and cache files.",
    include_common=["targets", "figr_config", "jobs", "verbose", "debug_log"],
    examples=[
        Example("./figr.sh clean", "Clean all build files"),
    ],
    key_options=[],
)

# HELP TOPICS

HELP_TOPICS = {
    "gpu": {
        "title": "GPU Configuration",
        "description": "How to configure GPU builds and runs",
    },
    "clusters": {
        "title": "Cluster Configuration",
        "description": "How to configure figr for different HPC clusters",
    },
    "batch": {
        "title": "Batch Job Submission",
        "description": "How to submit batch jobs with figr",
    },
    "debugging": {
        "title": "Debugging & Troubleshooting",
        "description": "Tips for debugging figr issues",
    },
}


# COMPLETE CLI SCHEMA

FIGR_CLI_SCHEMA = CLISchema(
    prog="./figr.sh",
    description="""\
Welcome to the figr master script. This tool automates and manages building, testing, \
running, and cleaning of figr in various configurations on all supported platforms. \
The README documents this tool and its various commands in more detail. To get \
started, run `./figr.sh build -h`.""",
    arguments=[
        Argument(
            name="help",
            short="h",
            help="Show help message",
            action=ArgAction.STORE_TRUE,
        ),
    ],
    commands=[
        BUILD_COMMAND,
        RUN_COMMAND,
        TEST_COMMAND,
        CLEAN_COMMAND,
    ],
    common_sets=[
        COMMON_TARGETS,
        COMMON_JOBS,
        COMMON_VERBOSE,
        COMMON_DEBUG_LOG,
        COMMON_GPUS,
        COMMON_FIGR_CONFIG,
    ],
    help_topics=HELP_TOPICS,
)


# DERIVED DATA (for use by other modules)

# Command aliases mapping (replaces COMMAND_ALIASES in user_guide.py)
COMMAND_ALIASES = {}
for cmd in FIGR_CLI_SCHEMA.commands:
    for alias in cmd.aliases:
        COMMAND_ALIASES[alias] = cmd.name


# Commands dict (replaces COMMANDS in user_guide.py)
def get_commands_dict():
    """Generate COMMANDS dict from schema for user_guide.py compatibility."""
    return {
        cmd.name: {
            "description": cmd.description or cmd.help,
            "alias": cmd.aliases[0] if cmd.aliases else None,
            "examples": [(e.command, e.description) for e in cmd.examples],
            "key_options": list(cmd.key_options),
        }
        for cmd in FIGR_CLI_SCHEMA.commands
    }


COMMANDS = get_commands_dict()
