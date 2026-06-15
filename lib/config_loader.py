import configparser
import os

class ResourceConfig:
    def __init__(self, config_file):
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Config file not found: {config_file}")
        self.cfg = configparser.ConfigParser()
        self.cfg.read(config_file)

    def get_params(self, step_name):
        """
        返回资源字典: {'cores': int, 'mem': int, 'partition': str}
        """
        if step_name not in self.cfg:
            # 如果找不到配置块，给予一个默认值，防止报错
            print(f"Warning: Section [{step_name}] not found in config, using defaults.")
            return {"cores": 1, "mem": 4096, "partition": None}
        
        return {
            "cores": self.cfg.getint(step_name, "threads"),
            "mem": self.cfg.getint(step_name, "memory_mb"),
            # 读取字符串，如果没填则返回 None (AnADAMA2 会使用集群默认分区)
            "partition": self.cfg.get(step_name, "partition", fallback=None)
        }