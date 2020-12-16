# CIP
Dataset and Code of CIP (Complementary-View Co-Interest Person Detection), published in ACM MM 2020.

```
@inproceedings{han2020cip,
  title={Complementary-View Co-Interest Person Detection}, 
  author={Han, Ruize and Zhao, Jiewen and Feng, Wei and Gan, Yiyang and Wan, Liang and Wang, Song},  
  year={2020},  
  booktitle={ACM International Conference on Multimedia}
}
```

## Introduction


<div align=center><img src="https://github.com/RuizeHan/CIP/blob/master/figs/example.png" width="450" height="360" alt="example"/><br/>
<div align= justify>
Fast and accurate identification of the co-interest persons, who draw joint interest of the surrounding people, plays an important role in social scene understanding and surveillance. Previous study mainly focuses on detecting co-interest persons from a single-view video. In this paper, we study a much more realistic and challeng-
ing problem, namely co-interest person (CIP) detection from multiple temporally-synchronized videos taken by the complementary and time-varying views. Specifically, we use a top-view camera, mounted on a flying drone at a high altitude to obtain a global view of the whole scene and all subjects on the ground, and multiple
horizontal-view cameras, worn by selected subjects, to obtain a local view of their nearby persons and environment details. We present an efficient top- and horizontal-view data fusion strategy to mapmultiplehorizontalviewsintotheglobaltopview.We then propose a spatial-temporal CIP potential energy function that jointly considers both intra-frame confidence and inter-frame consistency, thus leading to an effective Conditional Random Field (CRF) formulation. We also construct a complementary-view video dataset, which provides a benchmark for the study of multi-view co-interest person detection. Extensive experiments validate the effectiveness and superiority of the proposed method.

We show the CIP detection results on sample frames from the top-view and three horizontal-view videos in the following figure. Red and green boxes indicate the detected CIP and the ground truth, respectively. Frames with a solid blue star on the top-left corner indicate that no CIP is detected by our algorithm, e.g., they are drawn from the CIP’s egocentric video or the CIP is occluded in these frames. As shown in the figure, the proposed algorithm can detect CIP even if the CIP shows similar appearance and motion characteristics with other people. We can see that the proposed method can combine multiple horizontal views and a complementary top view for more comprehensive and accurate CIP detection.

![example](https://github.com/RuizeHan/CIP/blob/master/figs/cap_case.jpg)

Dataset: Link: https://pan.baidu.com/s/1BAM_MD1Lya5zf-O4OzJFnw  Password: CVIP 

To get the annotation, please contact the authors. The dataset is only used for academic research.

Authors: Ruize Han (han_ruize@tju.edu.cn); Jiewen Zhao (zhaojw@tju.edu.cn); Yiyang Gan (realgump@tju.edu.cn).
