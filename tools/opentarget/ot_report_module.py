"""Optional, conservative Russian-language Markdown interpretation of Open Targets features."""
from pathlib import Path
import pandas as pd

# All keys correspond to existing feature columns; no additional API calls.
EVIDENCE = {
    'Клинические данные': {
        'sm': [('approved_drug', 'Есть одобренный низкомолекулярный препарат, связанный с этой мишенью (Phase IV bucket).'), ('advanced_clinical', 'Низкомолекулярные препараты достигали фазы II или III.'), ('phase_1_clinical', 'Низкомолекулярные препараты достигали фазы I.')],
        'ab': [('approved_drug', 'Есть одобренный препарат на основе антител, направленный против мишени.'), ('advanced_clinical', 'Антительные препараты достигали фазы II или III.'), ('phase_1_clinical', 'Антительные препараты достигали фазы I.')],
        'pr': [('approved_drug', 'Есть одобренный препарат в модальности PROTAC согласно классификации платформы.'), ('advanced_clinical', 'PROTAC-модальность: клинический bucket фазы II/III.'), ('phase_1_clinical', 'PROTAC-модальность: клинический bucket фазы I.')],
        'oc': [('approved_drug', 'Есть одобренный препарат другой модальности.'), ('advanced_clinical', 'Другие модальности достигали фазы II или III.'), ('phase_1_clinical', 'Другие модальности достигали фазы I.')],
    },
    'Малые молекулы — экспериментальные и предсказанные предпосылки': {
        'sm': [
            ('structure_with_ligand', 'Есть экспериментальная структура комплекса мишени с малой молекулой (PDB).'),
            ('high_quality_ligand', 'Есть high-quality лиганды по критериям Open Targets/ChEMBL.'),
            ('high_quality_pocket', 'Предсказан высококачественный связывающий карман (DrugEBIlity ≥ 0.7).'),
            ('med_quality_pocket', 'Есть карман средней категории по DrugEBIlity (0–0.7).'),
            ('druggable_family', 'Белок принадлежит к семейству, которое считается потенциально druggable.'),
        ]
    },
    'PROTAC — предпосылки для направленной деградации': {
        'pr': [
            ('small_molecule_binder', 'Есть данные о связывающем маломолекулярном соединении из ChEMBL: важная предпосылка для дизайна деградера.'),
            ('literature', 'Мишень упоминается в вручную отобранной литературе о PROTAC; это не обязательно означает наличие работающего деградера.'),
            ('uniprot_ubiquitination', 'UniProt аннотирует убиквитинирование белка.'),
            ('database_ubiquitination', 'Сайты убиквитинирования зарегистрированы в специализированных базах.'),
            ('half_life_data', 'Измерения времени полужизни белка существуют, но сам период полужизни этот флаг не показывает.'),
        ]
    },
    'Антитела — предпосылки по локализации': {
        'ab': [
            ('uniprot_loc_high_conf', 'UniProt: высокоуверенное указание на доступную антителам локализацию.'),
            ('go_cc_high_conf', 'Gene Ontology: высокоуверенная аннотация доступной антителам локализации.'),
            ('human_protein_atlas_loc', 'Human Protein Atlas: данные о локализации в плазматической мембране.'),
            ('uniprot_loc_med_conf', 'UniProt: локализация со средней уверенностью.'),
            ('go_cc_med_conf', 'Gene Ontology: локализация со средней уверенностью.'),
            ('uniprot_sigp_or_tmhmm', 'Предсказаны сигнальный пептид или трансмембранные участки (SigP/TMHMM).'),
        ]
    }
}

ALL_KNOWN = {f'ot_{mod}_{suffix}' for group in EVIDENCE.values() for mod, pairs in group.items() for suffix, _ in pairs}


def is_true(x):
    if pd.isna(x):
        return False
    return str(x).strip().lower() in ('true', '1', '1.0')


def fmt_num(x):
    if pd.isna(x) or str(x).strip() == '':
        return 'нет данных'
    try:
        return f'{float(x):.3f}'
    except (ValueError, TypeError):
        return str(x)


def make_report(features_csv, output_file, max_genes=20):
    df = pd.read_csv(features_csv)
    if 'gene' not in df.columns:
        raise ValueError('CSV должен содержать столбец gene')
    n = df['gene'].dropna().nunique()
    if n > max_genes:
        raise ValueError(f'Отчёт отключён: {n} генов превышают лимит {max_genes}.')
    lines = [
        '# Open Targets — отчёт по терапевтической доступности мишеней', '',
        f'**Генов:** {n} · **Источник:** `{Path(features_csv).name}`', '',
        '> Отчёт интерпретирует только уже загруженные evidence buckets. `True` = удовлетворяет критерию категории, `False` = нет положительного evidence в этой категории, **а не** невозможность терапии.',
        '> Клинические стадии относятся к препаратам против мишени **при любом заболевании**, не обязательно JMML. Эти buckets отражают наличие исторических свидетельств, а не текущий статус исследования.',
        '> Разные buckets могут описывать пересекающиеся доказательства; их нельзя просто суммировать как независимые баллы.', '',
        '## Пояснения', '',
        '- **SM** — малые молекулы; **AB** — антитела; **PR** — PROTAC; **OC** — другие модальности.',
        '- **LoF O/E** — отношение наблюдаемого числа потерь функции к ожидаемому; **LOEUF** — верхняя граница интервала O/E. Низкие значения указывают на популяционную непереносимость LoF, но не количественно предсказывают токсичность препарата.',
        '- **Убиквитинирование** — косвенное свидетельство для PROTAC, **не подтверждение существования деградера**.',
        '- **Phase I, II/III, IV** — соответствующие категории Open Targets/ChEMBL; для названия препарата и показания нужен отдельный запрос known drugs.', '',
        '---', ''
    ]
    for _, row in df.iterrows():
        gene = str(row['gene'])
        eid = row.get('ensembl_id')
        status = str(row.get('ot_request_status', 'unknown'))
        lines += [f'## {gene}', '', f'**Ensembl:** {eid if pd.notna(eid) else "нет"} · **Статус получения данных:** {"успешно" if status == "ok" else "данные недоступны"}', '']
        for column, label in (('ot_approved_name', 'Название'), ('ot_pathways', 'Пути Reactome')):
            if pd.notna(row.get(column)):
                lines += [f'**{label}:** {row.get(column)}', '']
        if status != 'ok' or not is_true(row.get('ot_target_found')):
            lines += ['Нет подтверждённых данных из Open Targets: нельзя интерпретировать пропуски как отрицательные признаки.', '', '---', '']
            continue
        lines += [
            '### Genetic constraint', '',
            f'- LoF O/E: **{fmt_num(row.get("ot_lof_oe"))}**',
            f'- LOEUF: **{fmt_num(row.get("ot_lof_oe_upper"))}**',
            f'- Constraint score: **{fmt_num(row.get("ot_lof_constraint_score"))}** (сохраняется как предоставлен API; не вероятность токсичности)', '',
        ]
        for group_name, modalities in EVIDENCE.items():
            lines += [f'### {group_name}', '']
            any_positive = False
            for modality, pairs in modalities.items():
                present = [desc for suffix, desc in pairs if is_true(row.get(f'ot_{modality}_{suffix}'))]
                if not present:
                    continue
                any_positive = True
                if group_name == 'Клинические данные':
                    lines.append(f'**{modality.upper()}**')
                for desc in present:
                    lines.append(f'- {desc}')
                lines.append('')
            if not any_positive:
                lines += ['Положительных признаков в этой группе не найдено.', '']
        unknown = [col for col in df.columns if col.startswith(('ot_sm_', 'ot_ab_', 'ot_pr_', 'ot_oc_')) and col not in ALL_KNOWN and is_true(row.get(col))]
        if unknown:
            lines += ['### Дополнительные сведения', '', 'Есть дополнительные положительные категории, для которых пока не подготовлена проверенная текстовая интерпретация. Подробности остаются в CSV.', '']
        if pd.notna(eid):
            lines += [f'[Карточка Open Targets](https://platform.opentargets.org/target/{eid})', '']
        lines += ['---', '']
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text('\n'.join(lines), encoding='utf-8')
    return output_file