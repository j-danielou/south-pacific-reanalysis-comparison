# 🌊 SPRC (South-Pacific Reanalysis Comparison)

![Language](https://img.shields.io/badge/Language-R-blue)
![Status](https://img.shields.io/badge/Status-Active-success)

## 📖 À propos (About)

**SPRC** est une boîte à outils développée en **R** conçue pour faciliter la comparaison de produits de réanalyse climatique entre eux, ainsi qu'avec des données d'observation (in situ ou satellites). Le projet se concentre particulièrement sur la région du Pacifique Sud.

L'objectif de ce dépôt est de fournir des fonctions permettant d'évaluer la fiabilité et la précision des modèles de réanalyse à travers des statistiques robustes et des visualisations claires.

## ✨ Fonctionnalités (Features)

- **Métriques statistiques classiques :** Calcul automatisé de l'écart-type (SD), du coefficient de corrélation (r), de l'erreur quadratique moyenne (RMSE) et du biais.
- **Cartographie spatiale :** Génération de cartes pour visualiser facilement la distribution spatiale des variables et de leurs biais.
- **Diagrammes de Taylor :** Création de diagrammes de Taylor pour résumer graphiquement la correspondance entre les modèles et les observations (en combinant corrélation, variance et RMSE).

## 📂 Structure du dépôt

```text
south-pacific-reanalysis-comparison/
├── R/                     # Scripts R contenant les fonctions d'analyse et de visualisation
├── csv/                   # Répertoire contenant les fichiers de données (inputs/outputs)
├── plot_internship.Rproj  # Fichier de projet RStudio
└── README.md              # Documentation du projet
