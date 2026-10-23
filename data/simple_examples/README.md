# 小型示例

`raw/*.mps.gz` 是仓库原有 9 个手写 MPS 示例的无损 gzip 存档，涵盖界、
等式、混合不等式、自由变量、范围约束及退化换基等情况。

在仓库根目录执行（Python 3.9+，无需第三方包）：

```sh
python3 data/simple_examples/prepare.py --all
python3 data/simple_examples/prepare.py degenerate_pivot
```

默认在本目录生成同名 `.mps`，供现有测试直接读取；生成文件由 Git 忽略。
`--output-dir DIR` 可指定其他输出目录。输出通过临时文件生成，成功后替换。
`SHA256SUMS` 保存压缩原文件的校验值。
