import argparse
import json
import os


def load_config_file(parser: argparse.ArgumentParser, cfg_path: str) -> dict:
    if not os.path.isfile(cfg_path):
        parser.error(f"--config file does not exist: {cfg_path}")

    ext = os.path.splitext(cfg_path)[1].lower()
    try:
        with open(cfg_path, "r", encoding="utf-8") as handle:
            if ext == ".json":
                config = json.load(handle)
            elif ext in {".yml", ".yaml"}:
                try:
                    import yaml
                except ImportError:
                    parser.error("YAML config requested but PyYAML is not installed. Use JSON or install PyYAML.")
                config = yaml.safe_load(handle)
            else:
                parser.error(f"Unsupported config extension '{ext}'. Use .json, .yml, or .yaml.")
    except Exception as exc:
        parser.error(f"Failed to load config file {cfg_path}: {exc}")

    if config is None:
        config = {}
    if not isinstance(config, dict):
        parser.error("Config file must define a top-level object/dictionary of key-value pairs.")
    return config


def parse_args_with_config(
    parser: argparse.ArgumentParser,
    argv=None,
    defaults: dict | None = None,
    aliases: dict | None = None,
):
    cli_args = vars(parser.parse_args(argv))
    config_path = cli_args.pop("config", None)
    config_args = load_config_file(parser, config_path) if config_path else {}

    if aliases:
        for old_key, new_key in aliases.items():
            if old_key in config_args and new_key not in config_args:
                config_args[new_key] = config_args.pop(old_key)

    merged = {}
    if defaults:
        merged.update(defaults)
    merged.update(config_args)
    merged.update(cli_args)
    return argparse.Namespace(**merged)