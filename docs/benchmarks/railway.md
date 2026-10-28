# Railway 连续 LP 松弛 基准测试

简体中文 | [English](./railway_en.md) | [项目主页](../../README.md) | [基准测试汇总](./README.md)

由 `data/benchmark.py` 逐案例生成 MPS、执行求解器对比测试、清理 MPS 后更新。

## 测试与统计口径

- 计时为读取、求解和解导出的墙钟时间；不含 MPS 生成/删除、格式兼容处理、适配器公共初始化和独立审核。CLI 含进程启动，HiGHS 使用 Python binding，小案例开销有差异。
- PASS 要求求解器报告最优，原模型行约束、变量界和目标复算通过；有参考目标时再核对参考值（误差限为容差乘以 `1 + abs(reference)`）。无参考值显示 `—`，不代表独立最优性认证。
- 成功案例须全部重复为 PASS 且每次均在限时内；先取各案例重复耗时中位数，再计算普通几何平均（不移位、不对失败加惩罚）。空集合显示 `—`。
- 每个求解器分别与 lp_simplex 比较，仅使用双方共同成功的案例。几何平均耗时按“lp_simplex / 本求解器”列出；各行集合可能不同。相对耗时为 lp_simplex 耗时除以本求解器耗时，大于 1 表示本求解器更快；lp_simplex 行不列自比较统计。
- 成功数按案例计，超时与其他未通过数按重复运行计；统计仅属于对应运行的明确案例范围，进行中结果不构成全集结论。
- NOT_RUN / PENDING / INCOMPLETE 表示从未测试 / 本轮尚未测试 / 未完成全部重复。逐例单元格为“状态 / 秒数中位数 / 复算目标中位数”；非 PASS 的目标值仍须结合状态阅读。
- SKIPPED 表示按预设规则跳过，不计为超时或其他未通过，也不参与耗时统计。
- 同一机器、求解器、输入与参数不变时直接复用逐例测量；仅补测缺失或发生变化的项。对比统计按兼容的完整记录重新计算。

## 测试信息与统计

展示当前版本的逐例结果及对比统计。

<!-- benchmark-run:accepted-railway -->
| 测试信息 | 值 |
| --- | --- |
| 状态 / 进度 | 已完成；7/7 个案例 |
| 时间（UTC） | 2026-10-28T13:07:22.428932+00:00 |
| 对比求解器 | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| 案例范围 | railway 全集；rail507, rail516, rail582, rail2536, rail2586, rail4284, rail4872 |
| 系统 / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| 内存 / 逻辑 CPU | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex 版本 | b6f657753fe886c7de801ccff83be0ab4129a3ad；工作区有改动 |
| 构建 / 编译器 | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP 后端 | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| 算法（lp / scip / highs） | auto / default / choose |
| 重复 / 单次限时 / 容差 / 工作预算 | 1 / 3600 秒 / 1e-07 / 10000000 |
| 资源 | 串行；每个求解器单线程；未限制内存，不保证机器独占 |
| 原始记录 | 逐例测量与指纹保存在本报告中；临时文件已清理 |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK 算法 | dual simplex / primal simplex；presolve on |
| GLPK 跳过策略 | 任一案例超时后跳过后续测试；已跳过 3 次 |

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 7/7 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 7/7 | 0 | 0 | 7 | 34.2404 / 128.945 | 0.265543 |
| HiGHS | 7/7 | 0 | 0 | 7 | 34.2404 / 90.0054 | 0.380426 |
| Clp | 7/7 | 0 | 0 | 7 | 34.2404 / 53.5764 | 0.639095 |
| GLPK | 3/7 | 1 | 0 | 3 | 1.42014 / 3.19814 | 0.444052 |

复现本轮案例与参数（新建运行目录）：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset railway --cases rail507 rail516 rail582 rail2536 rail2586 rail4284 rail4872 --binary build/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 3600.0 --iterations 10000000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## 逐例最新结果

| 实例 | 规模（行 / 列 / 非零元） | 参考目标 | lp_simplex（状态 / 秒 / 目标） | SCIP / SoPlex（状态 / 秒 / 目标） | HiGHS（状态 / 秒 / 目标） | Clp（状态 / 秒 / 目标） | GLPK（状态 / 秒 / 目标） |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| rail507 | 507 / 63009 / 409349 | — | PASS / 1.626 / 172.145566677 | PASS / 5.551 / 172.145566677 | PASS / 4.188 / 172.145566677 | PASS / 3.866 / 172.145566677 | PASS / 3.994 / 172.145566677 |
| rail516 | 516 / 47311 / 314896 | — | PASS / 1.104 / 182 | PASS / 11.499 / 182 | PASS / 1.234 / 182 | PASS / 0.546 / 182 | PASS / 1.902 / 182 |
| rail582 | 582 / 55515 / 401708 | — | PASS / 1.596 / 209.712232957 | PASS / 7.596 / 209.712232957 | PASS / 3.052 / 209.712232957 | PASS / 2.594 / 209.712232957 | PASS / 4.306 / 209.712232957 |
| rail2536 | 2536 / 1081841 / 10993311 | — | PASS / 165.556 / 688.399223322 | PASS / 1104.127 / 688.399223322 | PASS / 1559.082 / 688.399223322 | PASS / 363.084 / 688.399223322 | TIMEOUT / — / — |
| rail2586 | 2586 / 920683 / 8008776 | — | PASS / 271.378 / 935.921820568 | PASS / 579.448 / 935.921820568 | PASS / 669.134 / 935.921820568 | PASS / 466.223 / 935.921820568 | SKIPPED / — / — |
| rail4284 | 4284 / 1092610 / 11279748 | — | PASS / 651.146 / 1054.05449127 | PASS / 1260.879 / 1054.05449127 | PASS / 1860.428 / 1054.05449127 | PASS / 1212.850 / 1054.05449127 | SKIPPED / — / — |
| rail4872 | 4872 / 968672 / 9244093 | — | PASS / 658.537 / 1509.63940162 | PASS / 1515.314 / 1509.63940162 | PASS / 1563.072 / 1509.63940162 | PASS / 1127.152 / 1509.63940162 | SKIPPED / — / — |

<!-- benchmark-cache:v1
[{"key":"40a140ab84b0203069accf2663aebab40a108ff2b858dfbe3b3677be43ef0d85","result":{"instance":"rail507","solver":"lp_simplex","repeat":1,"status":"PASS","wall_seconds":1.6261872090399265,"objective":172.14556667654884,"rows":507,"columns":63009,"nonzeros":409349,"audit_pass":true,"primal_residual":4.085620730620576e-14,"reference_error":null,"solver_status":"result:     state=0, code=0 (success)","cli_wall_seconds":1.6222071249503642,"cli_cpu_seconds":1.5478100000000001,"cli_user_seconds":1.517649,"cli_system_seconds":0.030161}},{"key":"79d80e0f98237cbd1d16c71fa9df7224e42d722bf1e9ec8f93c1daec7e4c4344","result":{"instance":"rail516","solver":"lp_simplex","repeat":1,"status":"PASS","wall_seconds":1.103668584022671,"objective":182.0,"rows":516,"columns":47311,"nonzeros":314896,"audit_pass":true,"primal_residual":2.220446049250313e-16,"reference_error":null,"solver_status":"result:     state=0, code=0 (success)","cli_wall_seconds":1.1001649159006774,"cli_cpu_seconds":1.091626,"cli_user_seconds":1.074431,"cli_system_seconds":0.017195}},{"key":"814f43e95e7a1fad719fe9fc2705ce576f2287f0bf9b8358bc1c63c8669141cd","result":{"instance":"rail582","solver":"lp_simplex","repeat":1,"status":"PASS","wall_seconds":1.5958261659834534,"objective":209.7122329571994,"rows":582,"columns":55515,"nonzeros":401708,"audit_pass":true,"primal_residual":8.193445921733655e-14,"reference_error":null,"solver_status":"result:     state=0, code=0 (success)","cli_wall_seconds":1.5907453331165016,"cli_cpu_seconds":1.564629,"cli_user_seconds":1.539175,"cli_system_seconds":0.025453999999999997}},{"key":"adcf4689f81e3bd60261d6400110f3a2e81fa317465f7ab8513577ea090d7364","result":{"instance":"rail2536","solver":"lp_simplex","repeat":1,"status":"PASS","wall_seconds":165.55621100007556,"objective":688.3992233222118,"rows":2536,"columns":1081841,"nonzeros":10993311,"audit_pass":true,"primal_residual":4.0782044385839e-16,"reference_error":null,"solver_status":"result:     state=0, code=0 (success)","cli_wall_seconds":165.50041391700506,"cli_cpu_seconds":158.89240700000002,"cli_user_seconds":157.24637,"cli_system_seconds":1.646037}},{"key":"27a8230ab4e0d6e677e5305760ea61bfd27be4411677182007c792a746610753","result":{"instance":"rail507","objective":172.145566677,"repeat":1,"solver":"scip","status":"PASS","wall_seconds":5.551}},{"key":"89d5208df31230dc5389b7e8ebae56672c59c66d30cb26cfa56f421334aa71f5","result":{"instance":"rail507","objective":172.145566677,"repeat":1,"solver":"highs","status":"PASS","wall_seconds":4.188}},{"key":"a7e2d6f18603d02e928148b2ba1fc611645201b3adcc5a86096539b86bbd7b59","result":{"instance":"rail507","objective":172.145566677,"repeat":1,"solver":"clp","status":"PASS","wall_seconds":3.866}},{"key":"84323826051a4b091ddabdc6b883bb8e47c5f7f8ea4933b886c2de01bbda5844","result":{"instance":"rail507","objective":172.145566677,"repeat":1,"solver":"glpk","status":"PASS","wall_seconds":3.994}},{"key":"f1dad37ce3bd741b533226d28394bd15c7fa5034f45075c67daae2803485deea","result":{"instance":"rail516","objective":182.0,"repeat":1,"solver":"scip","status":"PASS","wall_seconds":11.499}},{"key":"cb4a420e1f3417da6c03e3e9d3b67cf88be1f4ab029aecf506d87af9bad684e8","result":{"instance":"rail516","objective":182.0,"repeat":1,"solver":"highs","status":"PASS","wall_seconds":1.234}},{"key":"a89dec3b337f3d333be17f2a289f249825b5e257c0ee2c4355345d775c896c38","result":{"instance":"rail516","objective":182.0,"repeat":1,"solver":"clp","status":"PASS","wall_seconds":0.546}},{"key":"4a893d831bd29f3cc2ade06dad291feeb85a72e6dae3ce86f1e8602be1781e29","result":{"instance":"rail516","objective":182.0,"repeat":1,"solver":"glpk","status":"PASS","wall_seconds":1.902}},{"key":"ca4c0133587ec8ed5d0bb1a71f9c6d7ded2f8e99a8f896dd2a52e1ff58016268","result":{"instance":"rail582","objective":209.712232957,"repeat":1,"solver":"scip","status":"PASS","wall_seconds":7.596}},{"key":"515baa3579f64d32777d296d60b32e542673cbaf325a46168ae4d33e61607690","result":{"instance":"rail582","objective":209.712232957,"repeat":1,"solver":"highs","status":"PASS","wall_seconds":3.052}},{"key":"6ed93ccb581999efd5362cdade0fc37c4d983c84018e9efdbd949fa7d0a4091d","result":{"instance":"rail582","objective":209.712232957,"repeat":1,"solver":"clp","status":"PASS","wall_seconds":2.594}},{"key":"1d0488ca468bd248aa7484ce1dba2bb2010f288f355d97d7f6dde85bb20d7040","result":{"instance":"rail582","objective":209.712232957,"repeat":1,"solver":"glpk","status":"PASS","wall_seconds":4.306}},{"key":"fb82f3675d89bd99c5a40da46a0d5af0fe9ef5dc3e52b462aa4fd7a1e465e38d","result":{"instance":"rail2536","objective":688.399223322,"repeat":1,"solver":"scip","status":"PASS","wall_seconds":1104.127}},{"key":"5ca97ca97c86a73d433f21674d3dcdba3075feb984afeb1dc6d1a6298f6f4de8","result":{"instance":"rail2536","objective":688.399223322,"repeat":1,"solver":"highs","status":"PASS","wall_seconds":1559.082}},{"key":"533d5d35564268a7cb8a016f231da828509c46c4e66a90174811e217b6f937c5","result":{"instance":"rail2536","objective":688.399223322,"repeat":1,"solver":"clp","status":"PASS","wall_seconds":363.084}},{"key":"71cffd3bbb24cf8f1e24ac1bed150a28551d06ca0dc4d25a03d6178159f04da7","result":{"instance":"rail2536","repeat":1,"solver":"glpk","status":"TIMEOUT","wall_seconds":3600.0}},{"key":"368eb62e25e39bd27f33c83ea2698b9097eeea2070c12490a52fd4499fbce0cb","result":{"instance":"rail2586","objective":935.921820568,"repeat":1,"solver":"scip","status":"PASS","wall_seconds":579.448}},{"key":"17f14ab77720455f5c5e3916b0608ad4352d1e04345e5fd06cfcf4ac9ddc2422","result":{"instance":"rail2586","objective":935.921820568,"repeat":1,"solver":"highs","status":"PASS","wall_seconds":669.134}},{"key":"280d793ba27ff3a41c8f14f15eb70ae8a428171df9b267d96c1b68c6608f169c","result":{"instance":"rail2586","objective":935.921820568,"repeat":1,"solver":"clp","status":"PASS","wall_seconds":466.223}},{"key":"6c9726c42da4c3e7e47ef5df2a7445b825607d919bab981c8052fdba1d816fde","result":{"instance":"rail2586","repeat":1,"solver":"glpk","status":"SKIPPED","wall_seconds":null}},{"key":"1a9fdf83a04feaa03d40011a43cbb8f03d8c0f8723aafd17be5e16f259ac009f","result":{"instance":"rail4284","objective":1054.05449127,"repeat":1,"solver":"scip","status":"PASS","wall_seconds":1260.879}},{"key":"7cfdc4ba3623d3762eae3507db9cb36a58ec2d09ec4430fc79d642b50c421438","result":{"instance":"rail4284","objective":1054.05449127,"repeat":1,"solver":"highs","status":"PASS","wall_seconds":1860.428}},{"key":"ff9ec0933280ad6070411f87f5f993364e9acccab844e8d854d1d04f13c578b2","result":{"instance":"rail4284","objective":1054.05449127,"repeat":1,"solver":"clp","status":"PASS","wall_seconds":1212.85}},{"key":"27f5baca1a93d9a831e52f81aac12a635eea9703195093e56b517808299d4147","result":{"instance":"rail4284","repeat":1,"solver":"glpk","status":"SKIPPED","wall_seconds":null}},{"key":"581eee68debe2dbc17380d77c37126bd4e36bf7b76dbf8459fbca26b22118e02","result":{"instance":"rail4872","objective":1509.63940162,"repeat":1,"solver":"scip","status":"PASS","wall_seconds":1515.314}},{"key":"fe5899d2d999301f4dbe282ee90167c7480d5084df952f79b149dd41baac9cd1","result":{"instance":"rail4872","objective":1509.63940162,"repeat":1,"solver":"highs","status":"PASS","wall_seconds":1563.072}},{"key":"7e05dff3462fde153c9bb625633e502a617350fe1c28a84cc2a89d9282c14f83","result":{"instance":"rail4872","objective":1509.63940162,"repeat":1,"solver":"clp","status":"PASS","wall_seconds":1127.152}},{"key":"781d6848e7ce18deb55fa2194cf94c639ab168b4956da5e10d405b087783d49d","result":{"instance":"rail4872","repeat":1,"solver":"glpk","status":"SKIPPED","wall_seconds":null}},{"key":"7078b07748bcb23a71ff51b342037a142aa40e194912611353f322c41bb6fa46","result":{"instance":"rail2586","solver":"lp_simplex","repeat":1,"status":"PASS","wall_seconds":271.3779434999451,"objective":935.9218205677988,"rows":2586,"columns":920683,"nonzeros":8008776,"audit_pass":true,"primal_residual":3.738668580499663e-16,"reference_error":null,"solver_status":"result:     state=0, code=0 (success)","cli_wall_seconds":271.3311119999271,"cli_cpu_seconds":255.262245,"cli_user_seconds":253.064887,"cli_system_seconds":2.197358}},{"key":"220880bd55df7c2092acc1dbd0ec481aeea0aae6f575a8b2b4aaeba4409d0045","result":{"instance":"rail4284","solver":"lp_simplex","repeat":1,"status":"PASS","wall_seconds":651.1462121249642,"objective":1054.0544912748433,"rows":4284,"columns":1092610,"nonzeros":11279748,"audit_pass":true,"primal_residual":3.3306690738754696e-16,"reference_error":null,"solver_status":"result:     state=0, code=0 (success)","cli_wall_seconds":651.0517272078432,"cli_cpu_seconds":593.376694,"cli_user_seconds":587.431962,"cli_system_seconds":5.944732}},{"key":"aa90f20a04c52b718a24a5945ed7d487cc400b5c1e46fd33778f3c68bdbf9083","result":{"instance":"rail4872","solver":"lp_simplex","repeat":1,"status":"PASS","wall_seconds":658.5370465829037,"objective":1509.6394016217437,"rows":4872,"columns":968672,"nonzeros":9244093,"audit_pass":true,"primal_residual":3.1642466424841587e-12,"reference_error":null,"solver_status":"result:     state=0, code=0 (success)","cli_wall_seconds":658.4893603748642,"cli_cpu_seconds":632.969963,"cli_user_seconds":629.672844,"cli_system_seconds":3.297119}}]
-->
