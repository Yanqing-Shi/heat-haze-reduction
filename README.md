# Heat Haze

这是一个用 Visual Studio 创建的 C++ / OpenCV 小工具，用来手动修复照片里的局部热浪/空气扰动边缘。

## 项目里真正重要的文件

- `Heat Haze.sln`: Visual Studio 解决方案，双击它打开项目。
- `Heat Haze/main.cpp`: 程序入口，负责打开图片和窗口。
- `Heat Haze/WindowZoom.*`: OpenCV 窗口、缩放、拖动、鼠标交互、撤销。
- `Heat Haze/method.*`: 选区里的边缘检测和图像修复逻辑。
- `Heat Haze/PolyFit.*`: 对选中的边缘点做二次曲线拟合。
- `Heat Haze/DSC_1744.jpg`: 默认测试图片。

## 怎么运行

1. 用 Visual Studio 打开 `Heat Haze.sln`。
2. 选择 `x64` 平台。
3. 如果电脑里的 OpenCV 路径不是 `C:\Program Files\opencv\build`，在项目属性里改：
   - C/C++ -> 常规 -> 附加包含目录
   - 链接器 -> 常规 -> 附加库目录
   - 调试 -> 环境里的 `PATH`
4. 直接运行会打开 `DSC_1744.jpg`。也可以在命令行传入图片：

```powershell
.\x64\Debug\Heat Haze.exe "Heat Haze\DSC_2571.JPG"
```

## 在 VS Code 里运行

这个项目现在也带了 VS Code 配置：

- `Terminal -> Run Build Task`: 编译 Debug x64。
- `Run and Debug -> Debug Heat Haze`: 编译并启动调试。
- `Terminal -> Run Task -> Run Debug`: 编译并直接运行。

这些配置仍然调用 Visual Studio 的 MSBuild，所以电脑上需要安装 Visual Studio C++ 工具链和 OpenCV。默认路径是：

- Visual Studio: `C:\Program Files\Microsoft Visual Studio\2022\Community`
- OpenCV: `C:\Program Files\opencv\build`

如果你的安装路径不同，改 `.vscode/tasks.json`、`.vscode/launch.json` 和 `Heat Haze/Heat Haze.vcxproj` 里的路径。

## 操作

- 左键拖动：圈出要处理的热浪/边缘区域。
- 鼠标滚轮：缩放。
- 右键拖动：放大后平移。
- 鼠标中键：撤销上一次处理。

处理后的当前图像会写到 `result.jpg`。Canny、mask、poly 等中间调试图默认不再输出；需要看中间图时，把 `method.cpp` 里的 `kSaveDebugImages` 改成 `true`。

## 可以忽略的东西

- `.vs/`
- `x64/`
- `*.obj`, `*.pdb`, `*.tlog`

这些都是 Visual Studio 自动生成的缓存或编译结果，删掉后重新编译会再生成。

## 上传到 GitHub 前

这个文件夹应该作为独立仓库上传，不要再用 `C:\Users\todds` 作为 Git 仓库根目录。建议只上传源码、配置、README 和示例图片，不上传 `.vs/`、`x64/`、`*.pdb`、`*.obj` 等生成文件。
