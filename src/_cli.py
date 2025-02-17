import click
import importlib


class LazyGroup(click.Group):
    def __init__(self, *args, lazy_subcommands=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.lazy_subcommands = lazy_subcommands or {}

    def list_commands(self, ctx):
        base = super().list_commands(ctx)
        lazy = sorted(self.lazy_subcommands.keys())
        return base + lazy

    def get_command(self, ctx, cmd_name):
        if cmd_name in self.lazy_subcommands:
            return self._lazy_load(cmd_name)
        return super().get_command(ctx, cmd_name)

    def _lazy_load(self, cmd_name):
        import_path = self.lazy_subcommands[cmd_name]
        modname, cmd_object_name = import_path.rsplit(".", 1)
        mod = importlib.import_module(modname)
        cmd_object = getattr(mod, cmd_object_name)
        if not isinstance(cmd_object, click.BaseCommand):
            raise ValueError(
                f"Lazy loading of {import_path} failed by returning "
                "a non-command object"
            )
        return cmd_object


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "mri": "gmri2fem._cli.mri",
        "seg": "gmri2fem._cli.seg",
        "dti": "dti._cli.dti",
        "brainmeshing": "brainmeshing._cli.brainmeshing",
        "i2m": "i2m._cli.i2m",
        "dcm2nii-mixed": "gmri2fem._cli.dcm2nii_mixed",
    },
)
def gmri2fem():
    pass
