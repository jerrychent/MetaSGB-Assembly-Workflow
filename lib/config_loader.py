import configparser
import os


class ResourceConfig:
    """Load scheduler resources and tool options from resources.cfg."""

    # threads are used by tools; cores are requested from the scheduler.
    DEFAULT_RESOURCES = {
        "kneaddata": {"threads": 4, "cores": 4, "mem": 16000, "partition": None},
        "assembly": {"threads": 16, "cores": 16, "mem": 200000, "partition": None},
        "binning": {"threads": 16, "cores": 16, "mem": 32000, "partition": None},
        "drep": {"threads": 32, "cores": 32, "mem": 128000, "partition": None},
        "quantification": {"threads": 16, "cores": 16, "mem": 64000, "partition": None},
        "phylophlan_sgb": {"threads": 16, "cores": 16, "mem": 64000, "partition": None},
    }

    def __init__(self, config_file):
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Config file not found: {config_file}")

        # utf-8-sig accepts both BOM and non-BOM config files.
        self.cfg = configparser.ConfigParser(
            inline_comment_prefixes=("#", ";"),
            interpolation=None,
        )
        with open(config_file, encoding="utf-8-sig") as handle:
            self.cfg.read_file(handle)

    def _resolve_section(self, section, fallback_section=None):
        """Return the requested section or an optional fallback section."""
        if self.cfg.has_section(section):
            return section
        if fallback_section and self.cfg.has_section(fallback_section):
            return fallback_section
        return None

    def get_params(self, step_name, fallback_step=None):
        """Return tool threads and scheduler resources for one workflow step."""
        section = self._resolve_section(step_name, fallback_step)
        defaults = self.DEFAULT_RESOURCES.get(
            step_name,
            self.DEFAULT_RESOURCES.get(
                fallback_step,
                {"threads": 1, "cores": 1, "mem": 4096, "partition": None},
            ),
        )

        if section is None:
            print(f"Warning: Section [{step_name}] not found in config, using defaults.")
            return dict(defaults)

        threads = self.cfg.getint(section, "threads", fallback=defaults["threads"])
        cores = self.cfg.getint(section, "scheduler_cores", fallback=threads)
        partition = self.cfg.get(section, "partition", fallback=defaults["partition"])
        if isinstance(partition, str):
            partition = partition.strip() or None

        return {
            "threads": threads,
            "cores": cores,
            "mem": self.cfg.getint(section, "memory_mb", fallback=defaults["mem"]),
            "partition": partition,
        }

    def get(self, section, option, fallback=None):
        """Read a string option."""
        resolved = self._resolve_section(section)
        if resolved is None:
            return fallback
        value = self.cfg.get(resolved, option, fallback=fallback)
        return value.strip() if isinstance(value, str) else value

    def getint(self, section, option, fallback=None):
        """Read an integer option."""
        resolved = self._resolve_section(section)
        if resolved is None:
            return fallback
        return self.cfg.getint(resolved, option, fallback=fallback)

    def getfloat(self, section, option, fallback=None):
        """Read a float option."""
        resolved = self._resolve_section(section)
        if resolved is None:
            return fallback
        return self.cfg.getfloat(resolved, option, fallback=fallback)

    def getboolean(self, section, option, fallback=None):
        """Read a boolean option."""
        resolved = self._resolve_section(section)
        if resolved is None:
            return fallback
        return self.cfg.getboolean(resolved, option, fallback=fallback)
