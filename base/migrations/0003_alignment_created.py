# -*- coding: utf-8 -*-
# Generated by Django 1.9.7 on 2016-09-27 16:17
from __future__ import unicode_literals

import datetime
from django.db import migrations, models
from django.utils.timezone import utc


class Migration(migrations.Migration):

    dependencies = [
        ('base', '0002_seqrecord_display_order'),
    ]

    operations = [
        migrations.AddField(
            model_name='alignment',
            name='created',
            field=models.DateTimeField(
                auto_now_add=True, default=datetime.datetime(2016, 9, 27, 16, 17, 49, 116177, tzinfo=utc)
            ),
            preserve_default=False,
        ),
    ]
