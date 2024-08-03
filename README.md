# Bifrost-TKCM pack
This is the result of a personal study.  
Note : This pack contains only the Bifrost-2.10.0.0-Windows build.
  
これは個人的な学習のために制作したBifrostノードパックです。  
packにはWindows版のBifrost 2.10.0.0用ビルドだけが含まれます。

<img src="pack/compounds/icon/tkcm.png" width="80px">

# Install
Please set the following environment variables.  
環境変数の設定を行ってください
- BIFROST_LIB_CONFIG_FILES = <path_to_TKCM_packs>/pack/TKCM_config.json

# pack
下記のカスタムノードとサンプルグラフが含まれます。

## RBF solver
![rbf](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/afe7f340-fefc-4e8a-a475-86772174b10d)
https://vimeo.com/manage/videos/896561829

## InstantMeshes
![instantmeshes](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/e9edcbac-695c-4921-9cc8-b3083b4d6df3)

## Tags
- tag_mesh_border
- mesh_tag_promote

![292610317-0b3bb7a7-99ab-47dc-9096-ae0669d14f4e](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/86360edc-66e4-488d-baa1-f8e8da609078)

## Snap
![snap](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/435105a6-80fe-40d0-a112-7247b8d74331)

## Merge Mesh Point
![mergemeshpoint](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/cdb3043f-ba0a-4e16-b7fc-ab67d031f6ba)

## Quadrangulate Mesh
![quad](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/996773f4-a088-4874-bece-0b8cd58920a8)

## Dissolve
- dissolve_point
- dissolve_face
- dissolve_half_edge

![dp](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/622e3dec-2082-464c-a226-687212c78d8a)
![df](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/0f523079-013c-40a2-8fa4-50bd99a143be)
![de](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/ff2ecf91-a5e3-4db7-b4d8-48cbe81a46e3)

## Delete Unused Mesh Components
![dmc](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/b693f1f2-5dc5-4d86-b6df-9ed289c7c68e)

## Poly Edge Expand
![expand](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/8d698efa-b2ee-420c-ba81-1da7aa2a2455)

## Detach Mesh Components
![detach](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/d173f759-01ab-41e8-9045-0f0f0465eb0f)

## Split Edge
![splitedge](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/4d32d771-38b2-4177-a6fe-e61016ebf8b0)

## Convex Hull
- convex_hull
- VHACD 
![ch2](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/9392cbf0-1497-45a4-9bc6-ff31f1530448)
![ch](https://github.com/TKCM/Bifrost-TKCMPack/assets/13941074/f20bd294-057b-4322-9316-5348880180f1)

# Known issues
[split_edge]  
: Due to a bug in the transfer_properties_by_weighted_sum node, the transfer of point_component/face_vertex_component properties is not performed correctly.  
: transfer_properties_by_weighted_sum ノードの不具合により、point_component/face_vertex_component プロパティの転送が正しく実行されません。  

# pack develop environment  
packの開発環境  
- CMake 3.17
- VisualStudio 2019  
- BifrostSDK : Bifrost 2.10.0.0 Windows  
- Eigen : 3.4.0 (use only for RBF-solver build)

