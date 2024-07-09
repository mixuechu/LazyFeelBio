from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_cors import CORS
from sqlalchemy import Boolean, String, Integer

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://postgres:MXC921223mxc..@localhost/mr_db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['JSON_AS_ASCII'] = False

mr_db = SQLAlchemy(app)
CORS(app, supports_credentials=True)
app.app_context().push()


class MTask(mr_db.Model):
    __tablename__ = 'm_task'

    id = mr_db.Column(mr_db.String(255), primary_key=True)
    name = mr_db.Column(mr_db.String(255), nullable=True)
    state = mr_db.Column(mr_db.String(50), nullable=True)
    type = mr_db.Column(mr_db.String(50), nullable=True)
    input_params = mr_db.Column(mr_db.JSON, nullable=True)
    description = mr_db.Column(mr_db.Text, nullable=True)
    create_time = mr_db.Column(mr_db.DateTime, default=mr_db.func.current_timestamp())


class MGwasData(mr_db.Model):
    __tablename__ = 'm_gwas_data'
    gwas_id = mr_db.Column(mr_db.String(255), primary_key=True)
    name = mr_db.Column(mr_db.String(255), nullable=False)
    state = mr_db.Column(mr_db.String(50), default="MISSING")


class MRResult(mr_db.Model):
    __tablename__ = 'm_mr_result'
    id = mr_db.Column(Integer, primary_key=True, autoincrement=True)  # Primary key
    exposure_id = mr_db.Column(String(255), mr_db.ForeignKey('m_gwas_data.gwas_id'), nullable=False)
    outcome_id = mr_db.Column(String(255), mr_db.ForeignKey('m_gwas_data.gwas_id'), nullable=False)
    task_id = mr_db.Column(String(255), mr_db.ForeignKey('m_task.id'), nullable=False)
    finished = mr_db.Column(Boolean, default=False, nullable=False)

# class MEntity(mr_db.Model):
#     __tablename__ = 'm_entity'
#
#     id = mr_db.Column(mr_db.String(50), primary_key=True)
#     name = mr_db.Column(mr_db.String(255), nullable=False)
#     type = mr_db.Column(mr_db.String(50), nullable=False)
#     description = mr_db.Column(mr_db.Text)
#
#
# class MMrRelation(mr_db.Model):
#     __tablename__ = 'm_mr_relation'
#
#     id = mr_db.Column(mr_db.Integer, primary_key=True, autoincrement=True)
#     exposure_id = mr_db.Column(mr_db.String(50), mr_db.ForeignKey('m_gwas_data.gwas_id'), nullable=False)
#     outcome_id = mr_db.Column(mr_db.String(50), mr_db.ForeignKey('m_gwas_data.gwas_id'), nullable=False)
#     task_id = mr_db.Column(mr_db.String(50), mr_db.ForeignKey('m_task.id'), nullable=False)
#     effective_ids = mr_db.Column(mr_db.JSON)
#     method = mr_db.Column(mr_db.String(50))
#     p_value = mr_db.Column(mr_db.Float)
#     effect_size = mr_db.Column(mr_db.Float)
#     confidence_interval = mr_db.Column(mr_db.String)