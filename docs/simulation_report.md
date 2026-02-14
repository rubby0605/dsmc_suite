# 粒子模擬引擎的未來 — 從行星大氣到半導體製程

**塗翎 (Ruby Lin Tu)**
**2026-02-14**

---

## 為什麼現在是好時機

三件事正在同時發生：

1. **Europa Clipper 2030 年抵達木衛二** — NASA 史上最大行星探測任務，搭載 MASPEX 質譜儀和 UVS 紫外光譜儀，將首次直接量測 Europa 的水蒸氣大氣層。所有觀測數據的詮釋都需要大氣模型。

2. **半導體製程走向原子級** — 3nm 以下節點，蝕刻/鍍膜的 profile 控制需要粒子級模擬（PIC-DSMC），不再是連續流體就能搞定的。2025 年 GEC 會議上 GPU 加速 PIC-DSMC 是熱門題目。

3. **AI surrogate model 爆發** — 2024 年已有團隊用 neural network 取代行星大氣輻射傳輸模組，加速 147 倍，準確度 >99%。同樣的思路可以用在任何計算瓶頸上。

我手上有一套行星大氣 DSMC 模擬引擎（原碩博研究，已重構為 C + OpenMP），它的物理核心——粒子彈道、表面反應、氣相碰撞、MC 光線追蹤——跟半導體蝕刻/鍍膜模擬幾乎相同。

以下是幾個我認為值得做的方向。

---

## 方向 A：下一代行星大氣模擬

### 背景：領域正在轉變

行星大氣模擬正從「能跑就好」進入「要跑得快、跑得準、跑得多」的時代：

| 時期 | 代表工具 | 粒子數 | 特色 |
|------|---------|--------|------|
| 2000s | 個人 MATLAB code | ~10⁴ | 每個研究生自己寫 |
| 2010s | SPARTA (Sandia) | ~10⁶ | 開源 MPI 平行化 |
| 2020s | AMPS / PANTERA | ~10⁷ | GPU + PIC-DSMC |
| **未來** | **DSMC + ML** | **即時** | **surrogate model 取代暴力模擬** |

2025 年的關鍵進展：
- **Enceladus 羽流 DSMC 模擬**：用超算跑百萬粒子，發現噴發質量比舊估計少 20-40%
  （Mahieux et al., JGR Planets 2025）
- **行星大氣 surrogate model**：用 RNN 取代輻射傳輸模組，金星大氣模擬加速 147 倍
  （Teinturier et al., MNRAS 2024）
- **PANTERA**：Von Karman Institute 發布的開源 PIC-DSMC，Fortran + MPI

### 機會：我們能做什麼

**題目 1：Europa Clipper 觀測預測模型**

Europa Clipper 2030 年抵達木衛二，將進行 49 次飛越。
MASPEX 質譜儀會量測 Europa 極其稀薄的大氣成分。

但資料詮釋需要模型——探測器在不同飛越軌道上「應該」看到什麼密度？

→ 用 DSMC 跑 Europa 的 H₂O 大氣（改換物理常數即可），
   大規模參數掃描（噴發源位置、強度、表面溫度），
   建立 lookup table 或 ML surrogate model 供任務團隊使用。

**時程：** Europa Clipper 2030 才到，現在開始做模型正好趕上。
**發表：** 方法論文可先投，觀測比對論文等數據回來再投。

**題目 2：Multi-body 大氣比較研究**

同一套引擎，換參數就能跑不同天體：

| 天體 | 大氣成分 | 特殊物理 | 新數據來源 |
|------|---------|---------|-----------|
| Ceres | H₂O | 太陽活動觸發暫態大氣 | Dawn 遺產數據 |
| Europa | H₂O, O₂ | 木星磁層轟擊 | Clipper (2030) |
| Enceladus | H₂O | 南極噴泉 | Cassini 遺產 + 提案中的新任務 |
| Moon | Na, K, H₂O | 太陽風濺射 | 嫦娥、Artemis |
| Mercury | Na, Ca | 極端溫差 | BepiColombo (軌道中) |

→ 一篇比較論文就能涵蓋多個天體，引用數會比單一天體高。

**題目 3：ML Surrogate 取代 DSMC 計算瓶頸**

概念：不是取代整個模擬，而是取代最耗時的部分。

```
傳統：每個粒子跑完整軌道積分（10⁷ 時間步） → 慢

Surrogate：
  訓練集 = DSMC 跑 10000 個 (起始位置, 溫度, 速度) → (最終密度貢獻)
  推論   = neural network 秒級預測 → 快 1000 倍
  驗證   = 跟完整 DSMC 比對 → 確保精度
```

2024 年金星大氣的成功案例（147 倍加速）證明這條路可行。

---

## 方向 B：半導體蝕刻/鍍膜模擬

### 背景：為什麼粒子模擬在半導體越來越重要

隨著製程微縮到 3nm 以下，傳統連續流體模型失效：
- 蝕刻溝槽的 aspect ratio 超過 50:1
- 粒子平均自由程 > 特徵尺寸 → Knudsen number > 1
- 必須用粒子級模擬（DSMC / PIC-MCC）

2025 年的產業動態：
- **Lam Research** 公開宣布投資 digital twin 技術，用模擬取代實體原型
- **Tokyo Electron** 建立三層虛擬模型（製程→機台→整廠）
- **GEC 2025** 多篇 GPU 加速 PIC-DSMC 論文
- **AI-SPC** 在蝕刻製程中減少 40% 誤報、提升 1.7% 良率

### 物理對應：行星 ↔ 半導體

| 物理 | 行星大氣 | 半導體蝕刻/鍍膜 |
|------|---------|----------------|
| 粒子來源 | 表面昇華 | target 濺射 / plasma 解離 |
| 輸運 | 重力場彈道 | 無場/電場中飛行 |
| 碰撞 | 稀薄大氣 DSMC | 低壓腔體 DSMC |
| 表面反應 | sticking / 反彈 | 蝕刻 / 沉積 / 反射 |
| 光線追蹤 | 太陽遮擋 | 離子/中性粒子到溝槽底部的通量 |
| 溫度 | 表面能量平衡 | 基板溫度控制 |

核心物理引擎完全共用，只需要：
1. 改邊界條件（球面 → 平面/溝槽幾何）
2. 加電場模組（蝕刻用）
3. 加表面化學反應（蝕刻 vs 沉積概率）

### 機會：跟光電系合作

**現有開源工具的缺口：**

| 工具 | 特色 | 缺什麼 |
|------|------|--------|
| ViennaPS (Vienna) | Level-set + MC ray tracing，C++ | **沒有 DSMC 氣相碰撞** |
| SPARTA (Sandia) | 通用 DSMC，MPI | **不做表面演化** |
| PANTERA (VKI) | PIC-DSMC，Fortran | **偏航太，不做半導體** |
| 商業 TCAD | 功能完整 | **幾百萬台幣/年** |

→ **DSMC 碰撞引擎 + ViennaPS 的 level-set = 開源版完整蝕刻模擬器**
   這個組合目前不存在。

**具體合作方向：**

若廖瑩教授的實驗室有做以下任一項，就有合作空間：

| 實驗室製程 | 模擬能做的事 | 價值 |
|-----------|------------|------|
| PVD/Sputtering 鍍膜 | 預測膜厚均勻度 vs 靶距/氣壓 | 減少試片浪費 |
| PECVD 薄膜 | 模擬 step coverage | 優化溝槽填充 |
| RIE/ICP 蝕刻 | 預測蝕刻 profile 形狀 | 避免 undercut/bowing |
| ALD 原子層沉積 | 模擬前驅物輸運均勻度 | 改善大面積均勻性 |

**殺手級應用：AI + 模擬的 Digital Twin**

```
  實驗參數                              即時預測
(氣壓、功率、      ┌──────────────┐    (蝕刻 profile、
 溫度、氣體流量)  →│  ML Surrogate │→    膜厚分布、
                   │  (訓練自 DSMC) │     均勻度)
                   └──────────────┘
                          ↑
                    DSMC 模擬產生
                    數千組訓練資料
```

2025 年已有實例：
- ALD-GPR 模型用 MLP + Gaussian Process Regression 優化晶圓均勻度
- GlobalFoundries 用 AI 優化蝕刻製程，效率提升 5-10%
- NVIDIA 展示 AI-TCAD 加速，計算效率提升 3 個數量級

---

## 方向 C：方法論文 + 開源工具

不管最後做行星還是半導體，一篇方法論文都是好的起點：

> **"An Open-Source Parallel DSMC Framework for Planetary Exospheres
>   and Semiconductor Process Simulation"**

賣點：
1. 同一套物理引擎，兩個領域的應用
2. C + OpenMP，無外部依賴，容易移植
3. 跟 ViennaPS 互補（他們做表面，我們做氣相）
4. Benchmark 驗證（vs MATLAB 原版 + 解析解）
5. 開源 → 引用數高

投稿目標：
- Computer Physics Communications（偏方法）
- Journal of Computational Physics（偏數值）
- SoftwareX（偏開源工具）

---

## 我能帶什麼到桌上

| 能力 | 具體內容 |
|------|---------|
| 物理模擬 | 碩博班 DSMC 經驗 + 已完成 C 重構 |
| AI/ML | 業界 RD 經驗，可以做 surrogate model |
| 高效能計算 | OpenMP 已完成，可擴展 GPU/MPI |
| 軟體工程 | CMake、模組化架構、測試框架 |
| 時間 | Part-time side project，無畢業壓力 |

**對合作教授的好處：**
- 免費的模擬工具（vs 商業軟體）
- 模擬 + 實驗的論文比純實驗/純模擬更好發
- 跨領域合作本身就是賣點（行星科學方法 ↔ 半導體製程）
- 學生可以接著用、接著改

---

## 參考資料

### 行星科學
1. Europa Clipper Mission — NASA (2024 年 10 月發射，2030 抵達)
   https://science.nasa.gov/mission/europa-clipper/

2. MASPEX 質譜儀 — Europa 大氣量測核心儀器
   https://link.springer.com/article/10.1007/s11214-024-01061-6

3. Enceladus 羽流 DSMC 模擬 (Mahieux et al., JGR 2025)
   https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025JE009008

4. ML Surrogate 加速行星大氣模擬 147 倍 (MNRAS 2024)
   https://academic.oup.com/mnras/article/535/3/2210/7853137

5. 大氣外氣層模擬方法綜述 (Frontiers 2024)
   https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2024.1484360/full

6. SPARTA — Sandia 開源 DSMC (2025 更新)
   https://sparta.github.io/

7. PANTERA — VKI 開源 PIC-DSMC
   https://github.com/vonkarmaninstitute/pantera-pic-dsmc

8. Ceres 暫態大氣與太陽活動的關聯 — JPL
   https://www.jpl.nasa.gov/news/ceres-temporary-atmosphere-linked-to-solar-activity/

### 半導體製程模擬
9. ViennaPS — 開源半導體製程模擬框架 (2025 paper)
   https://www.sciencedirect.com/science/article/pii/S2352711025004194

10. GPU-Accelerated PIC-DSMC for Cu PVD (GEC 2025)
    https://archive.aps.org/gec/2025/fr5/3/

11. AI-TCAD 模擬加速 — NVIDIA (2025)
    https://developer.nvidia.com/blog/using-ai-physics-for-technology-computer-aided-design-simulations/

12. Digital Twin in Semiconductor — Lam Research
    https://newsroom.lamresearch.com/how-digital-twins-can-revolutionize-chipmaking

13. Digital Twin in Semiconductor — Tokyo Electron (2025)
    https://www.tel.com/blog/all/20250828_001.html

14. AI-SPC 蝕刻製程良率提升 (IJSRM 2025)
    https://ijsrm.net/index.php/ijsrm/article/view/6439

15. ML 優化 ALD 均勻度 (2025)
    https://www.sciencedirect.com/science/article/abs/pii/S2452414X25001025

16. Plasma Etching 的未來挑戰 (JVST-B 2024)
    https://pubs.aip.org/avs/jvb/article/42/4/041501/3297248
