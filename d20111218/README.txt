修士時の研究に用いた実験用のコードと書き捨てスクリプト

[実験用コード]
全結合 boltzmann machine
	FBM.cc
	FBM.hh

伊達准教の書いたサンプル手直し
	date_main.cc

倉田教授の論文より実装
	kurata_main.cc
	kurata2_main.cc
	kurata3_main.cc

結果データから実験中のアニメーションを再現するスクリプト
	fbm_playback.py
	fbm_playback2.py

option解析
	readcfg.py

実験結果のアニメーションを表示
	viewer.py

円の描画
	draw_object.py

エネルギー関数のdebug用
	energy.py

hintonダイアグラムにより素子の重みを出力
	hinton_dia1.py

ublas(線形代数ライブラリ)の出力結果を読み取れるように加工
	read_ublas.py

Makefileのターゲット
	date	伊達サンプルのコード生成
	kurata	倉田教授の論文より実装
	fbmview	アニメーションの再生


