from dataclasses import dataclass
from typing import Optional, List, Dict
from datetime import datetime

@dataclass
class MTask:
    """MTask"""
    id: int
    name: str
    state: str
    input_params: Dict
    description: str
    create_time: datetime
    type: str

@dataclass
class MGwasData:
    """MGwasData"""
    id: int
    type: str
    gwas_id: str
    name: str
    is_downloaded: bool
    file_path: Optional[str] = None

@dataclass
class MEntity:
    """MEntity"""
    id: int
    name: str
    type: str
    description: Optional[str] = None

@dataclass
class MMrRelation:
    """MMrRelation"""
    id: int
    exposure_id: int
    outcome_id: int
    effective_ids: List[int]
    method: Optional[str] = None
    p_value: Optional[float] = None
    effect_size: Optional[float] = None
    confidence_interval: Optional[str] = None
