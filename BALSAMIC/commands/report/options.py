"""Balsamic report specific options."""
import click

from BALSAMIC.constants.analysis import RULE_DELIVERY_MODES, RuleDeliveryMode
from BALSAMIC.constants.rules import DELIVERY_RULES

OPTION_RULES_TO_DELIVER = click.option(
    "-r",
    "--rules-to-deliver",
    multiple=True,
    type=click.Choice(DELIVERY_RULES),
    help="Specify the rules to deliver. The delivery mode selected via the --delivery-mode option.",
)

OPTION_DELIVERY_MODE = click.option(
    "-m",
    "--delivery-mode",
    type=click.Choice(RULE_DELIVERY_MODES),
    default=RuleDeliveryMode.APPEND.value,
    show_default=True,
    help=f"Append rules to deliver to the current delivery option ({RuleDeliveryMode.APPEND.value}) or deliver only "
    f"the ones specified ({RuleDeliveryMode.RESET.value})",
)

OPTION_PRINT_FILES = click.option(
    "-p",
    "--print-files",
    is_flag=True,
    default=False,
    show_default=True,
    help="Print list of analysis files. Otherwise only final count will be printed.",
)

OPTION_SHOW_ONLY_MISSING_FILES = click.option(
    "-m",
    "--show-only-missing",
    is_flag=True,
    default=False,
    show_default=True,
    help="Only show missing analysis files.",
)
