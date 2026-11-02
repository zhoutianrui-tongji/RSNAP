import ast
import pathlib
import unittest


ROOT = pathlib.Path(__file__).parents[1]


class StaticIntegrityTests(unittest.TestCase):
    def test_pipeline_modules_compile(self):
        for path in (ROOT / "script").glob("*.py"):
            with self.subTest(path=path.name):
                ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    def test_example_config_defines_all_pipeline_sections(self):
        tree = ast.parse((ROOT / "examples/config.example.py").read_text(encoding="utf-8"))
        names = {node.targets[0].id for node in tree.body if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)}
        self.assertTrue({"user_args", "qc_args", "at_args", "eer_args", "degi_args", "ti_args", "faa_args"} <= names)

    def test_readme_points_to_maintained_modules(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn("script/RSNAP.py", readme)
        self.assertIn("docs/architecture.md", readme)


if __name__ == "__main__":
    unittest.main()
