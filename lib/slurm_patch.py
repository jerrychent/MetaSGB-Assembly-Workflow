"""Slurm scheduler compatibility patches for anadama2."""


def patch_slurm_cpus_per_task():
    """Map anadama2 task cores to Slurm cpus-per-task.

    Anadama2 writes task cores as "#SBATCH -n", which requests multiple Slurm
    tasks. This workflow runs single-task, multithreaded tools, so the same
    value should be requested as "--cpus-per-task" instead.
    """
    from anadama2.grid import slurm

    if getattr(slurm.SLURMQueue, "_sgb_cpus_per_task_patch", False):
        return

    def append_optional_lines(template, value, prefix=""):
        """Append optional template lines from a string or iterable."""
        if not value:
            return
        if isinstance(value, str):
            template.append(prefix + value)
            return
        template.extend(prefix + item for item in value)

    def submit_template(self):
        template = [
            "#SBATCH -p ${partition}",
            "#SBATCH -N 1",
            "#SBATCH --ntasks=1",
            "#SBATCH --cpus-per-task=${cpus}",
            "#SBATCH -t ${time}",
            "#SBATCH --mem ${memory}",
            "#SBATCH -o ${output}",
            "#SBATCH -e ${error}",
        ]

        options = getattr(self, "options", None)
        image = getattr(self, "image", None)
        environment = getattr(self, "environment", None)

        append_optional_lines(template, options, "#SBATCH ")
        append_optional_lines(template, image, "singularity run ")
        append_optional_lines(template, environment)

        return template

    slurm.SLURMQueue.submit_template = submit_template
    slurm.SLURMQueue._sgb_cpus_per_task_patch = True
